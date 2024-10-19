import deal_data, json, os, argparse



if __name__ == "__main__":
    argparser = argparse.ArgumentParser()
    argparser.add_argument("--dataset", type=str, default="polybench",choices=["polybench","lore","tsvc"])
    args = argparser.parse_args()

    time = {}
    if args.dataset == "polybench":
        ori_dir = "./polybench"
        with open("./data/raw_data/rag_polybench_30_1010.json") as f:
            datas = json.load(f)
            for data in datas:
                filename = data["filename"]
                # print(filename)
                code = data["code"]
                file_path = os.path.join(ori_dir, filename + ".check.c")
                header_path = os.path.dirname(file_path)
                deal_data.polybench_compile(
                    file_path=file_path,
                    header_path=header_path,
                    polybench_header_path=os.path.join(ori_dir, "polybench/utilities"),
                    output_path=os.path.join(ori_dir, filename + ".checksum"),
                )

                data["ori_sum"], time[filename] = (
                    deal_data.polybench_checksum_time_test_exec(
                        os.path.join(ori_dir, filename + ".checksum")
                    )
                )
        with open("./data/raw_data/rag_polybench_30_1010.json", "w") as f, open(
            "./ori_result/polybench_ori_time.json", "w"
        ) as f1:
            json.dump(datas, f, indent=4)
            json.dump(time, f1)

    if args.dataset == "lore":
        ori_dir = "./lore/LORE_artificial"
        with open("./data/raw_data/rag_lore_50_1010.json") as f:
            datas = json.load(f)
        for data in datas:
            filename = data["filename"]
            print(filename)
            code = data["code"]
            file_path = os.path.join(ori_dir, filename + ".check.c")
            print(
                deal_data.lore_compile(
                    file_path=file_path,
                    lore_header_path=ori_dir,
                    output_path=os.path.join(ori_dir, filename + ".out"),
                )
            )
            data["ori_sum"], time[filename] = deal_data.lore_checksum_time_test_exec(
                os.path.join(ori_dir, filename + ".out")
            )
        with open("./data/raw_data/rag_lore_50_1010.json", "w") as f, open(
            "./ori_result/lore_ori_time.json", "w"
        ) as f1:

            json.dump(datas, f, indent=4)
            json.dump(time, f1)
    if args.dataset == "tsvc":
        ori_dir="./tsvc/cfiles"
        with open("./data/raw_data/rag_tsvc_85_1010.json")as f:

            for i,data in enumerate(datas):
                filename=data["filename"]

                print(filename)
                code=data["code"]
                deal_data.tsvc_compile(file_path=os.path.join(ori_dir,filename+".check.c"),output_path=os.path.join(ori_dir,filename+".checksum"))
                data["ori_sum"],time[filename]=deal_data.tsvc_checksum_time_test_exec(os.path.join(ori_dir,filename+".checksum"))
        with open("./data/raw_data/rag_tsvc_85_1010.json","w")as f,open("./ori_result/tsvc_ori_time.json",'w')as f1:
            json.dump(datas,f,indent=4)
            json.dump(time,f1)
