import json,argparse
import deal_data,os

if __name__=="__main__":
    result={}
    argparser=argparse.ArgumentParser()
    argparser.add_argument("--dataset",type=str,default="polybench",choices=["polybench","lore","tsvc"])
    args=argparser.parse_args()
    
    if args.dataset=="polybench":
        ori_dir="./polybench"

        result={}
        time={}
        def find_all_name(filename):
            with open("./polybench/polybench/utilities/benchmark_list","r")as f:
                lines=f.readlines()
                for line in lines:
                    if filename.split("_")[0] in line:
                        return os.path.join("polybench/polybench",line[2:-3])
                

        for root,dir_,files in os.walk("./data/pluto_code/polybench_pluto_code/pluto_code"):
            for file in files:
                if file.endswith(".c"):
                    filename=file.split(".")[0]
                    file_all_name=find_all_name(filename)
                    print(file_all_name)
                    print(deal_data.polybench_compile(file_path=os.path.join(root,file),
                                                            header_path=os.path.join(os.path.dirname(file_all_name)),
                                                            polybench_header_path=os.path.join(ori_dir,"utilities"),
                                                            output_path=os.path.join(root,filename+".time")))
                    result[file_all_name]=deal_data.polybench_checksum_time_test_exec(exec_path=os.path.join(root,filename+".time"),plutotime=True)
        with open ("./pluto_result/polybench_plutotime_result.json","w")as f:
            json.dump(result,f)  
            
    if args.dataset=="lore":
        result={}

        def find_all_name(filename):
            with open("./lore/benchmark_list","r")as f:
                lines=f.readlines()
                for line in lines:
                    if filename in line:
                        return line[:-1]
                

        for root,dir_,files in os.walk("./data/pluto_code/lore_pluto_code/pluto_code"):
            for file in files:
                if file.endswith(".c"):
                    
                    filename=file.split(".")[0]
                    if filename=='2_NPB_bt':
                        file_all_name='SCImark+NPB/2_NPB_bt'
                    elif filename=='2_NPB_bt1':
                        file_all_name='SCImark+NPB/2_NPB_bt1'
                    else:
                        file_all_name=find_all_name(filename)
                    print(filename)
                    print(file_all_name)
                    print(deal_data.lore_compile(file_path=os.path.join(root,file),lore_header_path="./lore/LORE_artificial",
                    output_path=os.path.join(root,file+".time")))
                    result[file_all_name]=deal_data.lore_checksum_time_test_exec(exec_path=os.path.join(root,file+".time"),plutotime=True)
        with open ("./pluto_result/lore_plutotime_result.json","w")as f:
            json.dump(result,f)  
            
    if args.dataset=="tsvc":
        for root,dir_,files in os.walk("./data/pluto_code/tsvc_pluto_code/pluto_code"):
            for file in files:
                if file.endswith(".c"):
                    print(file)
                    print(deal_data.tsvc_compile(file_path=os.path.join(root,file),
                                                            output_path=os.path.join(root,file+".time")))
                    result[file]=deal_data.tsvc_checksum_time_test_exec(exec_path=os.path.join(root,file+".time"),plutotime=True)
        with open ("./pluto_result/tsvc_plutotime_result.json","w")as f:
            json.dump(result,f)  