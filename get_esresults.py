import argparse, insert_esdata, re
from deal_data.llm_response import llm_api
import deal_data
import json, pandas as pd




if __name__ == "__main__":
    argparser = argparse.ArgumentParser()
    argparser.add_argument("--dataset", type=str, default="polybench",choices=["polybench","lore","tsvc"])
    argparser.add_argument("--es_host", type=str, default="http://10.214.241.234:9210")
    argparser.add_argument("--index_name", type=str, default="")
    argparser.add_argument("--file_path", type=str, default="./mapping7.json")
    data_path="./data"

    args = argparser.parse_args()
    dataset = args.dataset
    if dataset == "tsvc":
        file_path = "./data/raw_data/rag_tsvc_85_1010.json"
        save_path = "./es_result/tsvc_es.xlsx"
    elif dataset == "polybench":
        file_path = "./data/raw_data/rag_polybench_30_1010.json"
        save_path = "./es_result/polybench_es.xlsx"
    elif dataset == "lore":
        file_path = "./data/raw_data/rag_lore_50_1010.json"
        save_path = "./es_result/lore_es.xlsx"
    else:
        file_path = "./data/raw_data/sampled_output.json"
        save_path = "./es_result/poly_code_es.xlsx"
    result_list = []
    with open(file_path, "r") as searchcodes:
        searchcodes = json.load(searchcodes)
        for i, data in enumerate(searchcodes):
            search_code, basename = data["code"], data["filename"]
            query_code=data
            query_code.pop("ori_sum")
            results = deal_data.use_rag(query_code)
            # 打印搜索结果
            print("search_code:\n",search_code)
            for hit in results["hits"]["hits"]:
                print(f"Code:\n {hit['_source']['code']}")
                print(f"Optimized Code:\n {hit['_source']['opt_code']}")
                # print(f"Loop Levels: {hit['_source']['loop_levels']}")
                print(f"Score: {hit['_score']}")
                print("-----")
                result_list.append(
                    {
                        "score": hit["_score"],
                        "_source": hit["_source"],
                        "filename": basename,
                    }
                )
    data = pd.DataFrame(result_list)
    data.to_excel(save_path, index=False)
