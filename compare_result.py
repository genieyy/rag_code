import json, argparse, deal_data
from collections import defaultdict

datas = defaultdict(lambda: defaultdict(int))

if __name__ == "__main__":
    argparser = argparse.ArgumentParser()
    argparser.add_argument("--dataset", type=str, default="polybench",choices=["polybench","lore","tsvc"])
    args=argparser.parse_args()

    datas = defaultdict(lambda: defaultdict(int))
    deal_data.extract_result(args.dataset)
    
    #############save all_results############
    with open(
        f"./{args.dataset}/results_rollback_timefeedback_all/rag_all_results.json",
        "r",
    ) as f:
        all_results = json.load(f)
        for file in all_results["compile1"]:
            values1 = list(all_results["compile1"][file].values())
            for i in range(len(values1) - 1):
                datas["compile1"][str(i + 1)] += values1[i]
            datas["compile1"]["final"] += values1[-1]
        for file in all_results["compile2"]:
            values2 = list(all_results["compile2"][file].values())
            for i in range(len(values2) - 1):
                datas["compile2"][str(i + 1)] += values1[i]
            datas["compile2"]["final"] += values2[-1]
        for file in all_results["applied_passes"]:
            values3 = list(all_results["applied_passes"][file].values())
            for i in range(len(values3) - 1):
                datas["applied_passes"][str(i + 1)] += values3[i]
            datas["applied_passes"]["final"] += values3[-1]
        for file in all_results["checksum_passes"]:
            values4 = list(all_results["checksum_passes"][file].values())
            for i in range(len(values4) - 1):
                datas["checksum_passes"][str(i + 1)] += values4[i]
            datas["checksum_passes"]["final"] += values4[-1]
        for file in all_results["elemwise_passes"]:
            values5 = list(all_results["elemwise_passes"][file].values())
            for i in range(len(values5) - 1):
                datas["elemwise_passes"][str(i + 1)] += values5[i]
            datas["elemwise_passes"]["final"] += values5[-1]
    with open(f"./{args.dataset}_compare_other.json", "w") as f:
        json.dump(datas, f)
    
    with open(
        f"./{args.dataset}/results_rollback_timefeedback_all/norag_all_results.json",
        "r",
    ) as f:
        all_results = json.load(f)
        for file in all_results["compile1"]:
            values1 = list(all_results["compile1"][file].values())
            for i in range(len(values1) - 1):
                datas["compile1"][str(i + 1)] += values1[i]
            datas["compile1"]["final"] += values1[-1]
        for file in all_results["compile2"]:
            values2 = list(all_results["compile2"][file].values())
            for i in range(len(values2) - 1):
                datas["compile2"][str(i + 1)] += values1[i]
            datas["compile2"]["final"] += values2[-1]
        for file in all_results["applied_passes"]:
            values3 = list(all_results["applied_passes"][file].values())
            for i in range(len(values3) - 1):
                datas["applied_passes"][str(i + 1)] += values3[i]
            datas["applied_passes"]["final"] += values3[-1]
        for file in all_results["checksum_passes"]:
            values4 = list(all_results["checksum_passes"][file].values())
            for i in range(len(values4) - 1):
                datas["checksum_passes"][str(i + 1)] += values4[i]
            datas["checksum_passes"]["final"] += values4[-1]
        for file in all_results["elemwise_passes"]:
            values5 = list(all_results["elemwise_passes"][file].values())
            for i in range(len(values5) - 1):
                datas["elemwise_passes"][str(i + 1)] += values5[i]
            datas["elemwise_passes"]["final"] += values5[-1]
    with open(f"./{args.dataset}_norag_compare_other.json", "w") as f:
        json.dump(datas, f)
        
    ###########################count number that final is the best among all methods####
    with open(
        f"./{args.dataset}/results_rollback_timefeedback_all/rag_all_results.json",
        "r",
    ) as f:
        all_results = json.load(f)
        cnt = 0
        for file in all_results["run_times"]:
            if (
                min(list(all_results["run_times"][file].values())[:-1])
                > list(all_results["run_times"][file].values())[-1]
            ):
                cnt += 1
        print("number that final is the best among all methods:",cnt)  

    ###############compare rag with norag####################################################
    datas = {}
    with open(
        f"./{args.dataset}/results_rollback_timefeedback_all/rag_all_results.json",
        "r",
    ) as f:
        all_results = json.load(f)

        for file in all_results["run_times"]:
            datas[file] = min(list(all_results["run_times"][file].values()))
            
    with open(f"./{args.dataset}_comparetime.json", "w") as f:
        f.write(json.dumps(datas))

    datas = {}
    with open(
        f"./{args.dataset}/results_rollback_timefeedback_all/norag_all_results.json",
        "r",
    ) as f:
        all_results = json.load(f)

        for file in all_results["run_times"]:
            datas[file] = min(list(all_results["run_times"][file].values()))
    with open(f"./{args.dataset}_comparetime_norag.json", "w") as f:
        f.write(json.dumps(datas))

    cnt=0
    with open (
        f"./{args.dataset}_norag_comparetime.json",
        "r",
    ) as f,open(
        f"./{args.dataset}_comparetime.json",
        "r",
    ) as f1:
        norag=json.load(f)
        rag=json.load(f1)
        for i in rag:
            if rag[i]<norag[i] and rag[i]!=float("inf") and norag[i]!=float("inf"):
                cnt+=1
                print(i)
                print(rag[i],norag[i])
        print("number that llm with rag overcomes norag:",cnt)
            
            
        
        #################compare with pluto#############################
    if args.dataset == "polybench":
        with open("./polybench_comparetime.json", "r") as f1, open(
            "./pluto_result/polybench_plutotime_result.json", "r"
        ) as f2:
            times = dict(json.load(f1))
            cnt = 0
            plutotime = dict(json.load(f2))
            for i, file in enumerate(times):
                if "polybench/"+file in plutotime:
                    if times[file] < plutotime["polybench/"+file]:
                        print(file)
                        cnt += 1
            print("number that llm with rag overcomes pluto:",cnt)
        ###########################draw_compare############################
        # with open("./polybench_comparetime_norag.json", "r")as noragtime,open("./polybench_plutotime_result.json")as plutotime,open("./polybench_comparetime.json","r")as ragtime,open("./polybench_ori_time.json")as oritime:
        #     pt=json.load(plutotime)
        #     nt=json.load(noragtime)
        #     rt=json.load(ragtime)
        #     ori=json.load(oritime)
        #     p_o=[]
        #     n_o=[]
        #     r_o=[]
        #     for file in pt:
        #         file_=file[10:]
        #         p_o.append(pt[file]/ori[file_])
        #         r_o.append(rt[file_]/ori[file_])
        #         n_o.append(nt[file_]/ori[file_])
            
        #     print(p_o,"\n",r_o,"\n",n_o)
        #     print(list(rt.keys()))
        

    elif args.dataset == "lore":
        with open("./lore_comparetime.json","r")as f1,open("./pluto_result/lore_plutotime_result.json","r")as f2:
            times=dict(json.load(f1))
            cnt=0
            plutotime=dict(json.load(f2))
            for i,file in enumerate(times):
                
                if file in plutotime:
                    if times[file]<plutotime[file]:
                        cnt+=1
                        print(file)
            print("number that llm with rag overcomes pluto:",cnt)
    else:
        with open("./tsvc_comparetime.json","r")as f1,open("./pluto_result/tsvc_plutotime_result.json","r")as f2:
            times=dict(json.load(f1))
            cnt=0            
            plutotime=dict(json.load(f2))
            for i,file in enumerate(times):               
                if file+".pluto.c" in plutotime:
                    if times[file]<plutotime[file+".pluto.c"]:
                        cnt+=1
            print("number that llm with rag overcomes pluto:",cnt)
                
