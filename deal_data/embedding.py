import torch, json, os
import torch.nn.functional as F
import torch.nn as nn
from torch import Tensor
from transformers import AutoTokenizer, AutoModel
import torch.distributed as dist
import torch.multiprocessing as mp
from torch.cuda.amp import autocast


def last_token_pool(last_hidden_states: Tensor, attention_mask: Tensor) -> Tensor:
    left_padding = attention_mask[:, -1].sum() == attention_mask.shape[0]
    if left_padding:
        return last_hidden_states[:, -1]
    else:
        sequence_lengths = attention_mask.sum(dim=1) - 1
        batch_size = last_hidden_states.shape[0]
        return last_hidden_states[
            torch.arange(batch_size, device=last_hidden_states.device), sequence_lengths
        ]


def batch_embedding(code_list):
    tokenizer = AutoTokenizer.from_pretrained(
        "/home/cyy/models/Salesforce/SFR-Embedding-Mistral"
    )
    model = AutoModel.from_pretrained(
        "/home/cyy/models/Salesforce/SFR-Embedding-Mistral"
    ).half()
    model = nn.DataParallel(model).to("cuda")

    max_length = 1024
    all_embeddings = torch.empty((1, 4096))
    batch_size = 196
    # Split code_list into batches
    for i in range(0, len(code_list), batch_size):
        batch_code_list = code_list[i : i + batch_size]

        batch_dict = tokenizer(
            batch_code_list,
            max_length=max_length,
            padding=True,
            truncation=True,
            return_tensors="pt",
        ).to("cuda")
        with autocast():
            with torch.no_grad():
                outputs = model(**batch_dict)
            embeddings = last_token_pool(
                outputs.last_hidden_state, batch_dict["attention_mask"]
            ).to("cpu")

        all_embeddings = torch.concat((all_embeddings, embeddings), dim=0)

        # 清理显存
        del batch_dict, outputs, embeddings
        torch.cuda.empty_cache()

    return all_embeddings.tolist()[1:]

    # world_size = 4
    # manager = mp.Manager()
    # return_dict = manager.dict()

    # mp.spawn(batch_embedding_rank, args=(world_size, code_list, return_dict), nprocs=world_size, join=True)

    # if 'embeddings' in return_dict:
    #     embeddings = return_dict['embeddings']
    #     return embeddings


def batch_embedding_rank(rank, world_size, code_list, return_dict):
    batch_size = 2  # 可以增加批处理大小
    os.environ["MASTER_ADDR"] = "localhost"
    os.environ["MASTER_PORT"] = "12355"

    dist.init_process_group("nccl", rank=rank, world_size=world_size)

    tokenizer = AutoTokenizer.from_pretrained(
        "/home/cyy/models/Salesforce/SFR-Embedding-Mistral"
    )
    model = (
        AutoModel.from_pretrained(
            "/home/cyy/models/Salesforce/SFR-Embedding-Mistral"
        )
        .half()
        .to(rank)
    )
    model = nn.parallel.DistributedDataParallel(model, device_ids=[rank])

    print(
        f"device {rank}: {torch.cuda.memory_summary(torch.cuda.current_device(), abbreviated=False)}"
    )

    max_length = 1024
    all_embeddings = []

    # Split code_list into batches
    for i in range(0, len(code_list), batch_size):
        batch_code_list = code_list[i : i + batch_size]

        batch_dict = tokenizer(
            batch_code_list,
            max_length=max_length,
            padding=True,
            truncation=True,
            return_tensors="pt",
        ).to(rank)

        print(
            f"device {rank}: {torch.cuda.memory_summary(torch.cuda.current_device(), abbreviated=False)}"
        )

        with autocast():
            outputs = model(**batch_dict)
            embeddings = last_token_pool(
                outputs.last_hidden_state, batch_dict["attention_mask"]
            ).to("cpu")

        all_embeddings.append(embeddings)

        # 清理显存
        del batch_dict, outputs, embeddings
        torch.cuda.empty_cache()

    all_embeddings = torch.cat(all_embeddings, dim=0)

    gathered_embeddings = [torch.zeros_like(all_embeddings) for _ in range(world_size)]
    dist.all_gather(gathered_embeddings, all_embeddings)

    gathered_embeddings = torch.cat(gathered_embeddings, dim=0)

    if rank == 0:
        return_dict["embeddings"] = gathered_embeddings

    dist.destroy_process_group()


def get_embedding(code):

    return batch_embedding([code])[0]


if __name__ == "__main__":
    with open("./sampled_output.json", "r") as f1, open(
        "./sampled_output_embedding.json", "w"
    ) as f2:
        data = json.load(f1)
        code_lst = []
        opt_code_lst = []
        for i in data:
            code_lst.append(i["code"])
            opt_code_lst.append(i["opt_code"])

        code_embed_lst = batch_embedding(code_lst)
        opt_code_embed_lst = batch_embedding(opt_code_lst)

        data_embedding = []
        for i in data:
            data_embedding.append(i)
            data_embedding[-1]["embedding"] = code_embed_lst.pop(0)
            data_embedding[-1]["opt_embedding"] = opt_code_embed_lst.pop(0)
        json.dump(data_embedding, f2)

    with open("./remaining_output.json", "r") as f1, open(
        "./remaining_output_embedding.json", "w"
    ) as f2:
        data = json.load(f1)
        code_lst = []
        opt_code_lst = []
        for i in data:
            code_lst.append(i["code"])
            opt_code_lst.append(i["opt_code"])
        code_embed_lst = batch_embedding(code_lst)
        opt_code_embed_lst = batch_embedding(opt_code_lst)

        data_embedding = []
        for i in data:
            data_embedding.append(i)
            data_embedding[-1]["embedding"] = code_embed_lst.pop(0)
            data_embedding[-1]["opt_embedding"] = opt_code_embed_lst.pop(0)
        json.dump(data_embedding, f2)
