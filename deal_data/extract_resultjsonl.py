import os
import json
from collections import defaultdict

def extract_result(dataset):
    all_results=defaultdict(lambda:defaultdict(int))
    if dataset not in ["polybench","lore", "tsvc"]:
        raise ValueError("dataset not supported")

    with open(f"./{dataset}/results_rollback_timefeedback_all/multigenerate_result.jsonl", "r") as f:
        for line in f:
            data = json.loads(line)
            if data["name"]=="root":
                message=data["message"]
                filename=list(message["c1"].keys())[0][:-2]
                all_results["compile1"][filename]=data["message"]["c1"]
                all_results["compile2"][filename]=data["message"]["c2"]
                all_results["applied_passes"][filename]=data["message"]["a"]
                all_results["checksum_passes"][filename]=data["message"]["check"]
                all_results["elemwise_passes"][filename]=data["message"]["elemcheck"]
                all_results["run_times"][filename]=data["message"]["run"]
            
    with open(f"./{dataset}/results_rollback_timefeedback_all/rag_all_results.json","w")as f:
        json.dump(all_results,f)
        
    all_results=defaultdict(lambda:defaultdict(int))

    with open(f"./{dataset}/results_rollback_timefeedback_all/multigenerate_result_norag.jsonl", "r") as f:
        for line in f:
            data = json.loads(line)
            if data["name"]=="root":
                message=data["message"]
                filename=list(message["c1"].keys())[0][:-2]
                all_results["compile1"][filename]=data["message"]["c1"]
                all_results["compile2"][filename]=data["message"]["c2"]
                all_results["applied_passes"][filename]=data["message"]["a"]
                all_results["checksum_passes"][filename]=data["message"]["check"]
                all_results["elemwise_passes"][filename]=data["message"]["elemcheck"]
                all_results["run_times"][filename]=data["message"]["run"]
            
    with open(f"./{dataset}/results_rollback_timefeedback_all/norag_all_results.json","w")as f:
        json.dump(all_results,f)
if __name__ == "__main__":
    extract_result("tsvc")
            
            