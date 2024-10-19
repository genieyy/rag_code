# -*- coding: utf-8 -*-
import argparse, insert_esdata, re, os
from deal_data.llm_response import llm_api
import deal_data
import json
from collections import defaultdict
import log, pathlib
from typing import Union
from deal_data.data_deal import extract_and_concatenate,extract_code_and_text



def main(args):
    

    ori_dir = "./lore/LORE_artificial"


    results_dir = "./lore/results_rollback_timefeedback_all"
    os.makedirs(results_dir, exist_ok=True)

    if args.use_rag:
        lore_record_path = "./lore/results_rollback_timefeedback_all/multigenerate_result.jsonl"
        opt_dir="./lore/LORE_artificial"
    else:
        lore_record_path = "./lore/results_rollback_timefeedback_all/multigenerate_result_norag.jsonl"
        opt_dir="./lore/norag/LORE_artificial"
        os.makedirs(opt_dir, exist_ok=True)
    logger = log.setup_logging(record_path=lore_record_path)  # get a logger object
    #######################some util functions#######

    def first_prompt(search_code: str, use_rag: bool, results=None):
        """
        first throw source_code into llm and get an opt_code from it.
        if use_rag is true, use rag to search and prompt making prompt from it.
        else just prompt making prompt from it.
        """
        prompt=""
        if use_rag:
            for i, hit in enumerate(results):
                if hit["_score"] < args.threshold:
                    continue
                prompt += f"""
                
    Example:

    // original code
    {hit["_source"]["code"]}

    // optimized code
    {hit["_source"]["opt_code"][:1000]}"""
            prompt+=f"""

    Please analyze what meaning-preserving loop transformation methods are used in above examples, and tell me what you learn.

    please use appropriate methods you learn from examples to improve its performance:
    {search_code}

    1. Provide one optimized code.
    2. Do not include the original C program in your response.
    3. Do not define new function. 
    4. Existing variables do not need to be redefined. If you generate a new variable, please use the double type.
    5. Put your code in the code block(like ```c for i=1 .. ```).

        """


        else:
            prompt = f"""as a compiler, given the C program below, improve its performance using meaning-preserving transformation:
    {search_code}

    1. Provide one optimized code.
    2. Do not include the original C program in your response.
    3. Do not define new function. 
    4. Existing variables do not need to be redefined. If you generate a new variable, please use the double type.
    5. Put your code in the code block(like ```c for i=1 .. ```).
    """

        return prompt

    ####################################################
    def get_newoptcode_from_llm(messages: list, error: str, last_code, again: str = ""):
        """
        throw error into llm and regenerate a new optcode
        """

        massage = f"""this optimized version:\n{last_code}\n did a wrong transformation from the source code{ again}, resulting in a compilation error. 
    This is the compiler error message:{error}
    Please check the optimized version and regenerate it.
    Do not define a new function. The original variable does not need to be redefined. If you generate a new variable, please use the double type.
    put your code in the markdown code block."""
        messages = messages.copy()
        messages.append({"role": "assistant", "content": last_code})
        messages.append({"role": "user", "content": massage})
        new_response = llm_api[args.llm](messages)
        try:
            new_opt_code = extract_and_concatenate(new_response)
            print("new_opt_code: \n", new_opt_code)
            return new_opt_code
        except:
            return new_response


    def get_final_code_from_llm(search_code: str, feedback_list: list) -> str:
        """
        final code after time feedback to llm
        """
        if not feedback_list:  #########和no_rag同一个prompt
            prompt = f"""as a compiler, given the C program below,improve its performance using meaning-preserving transformation.
    Do not define a new function. The original variable does not need to be redefined. If you generate a new variable, please use the double type.Put your code in the code block(like ```c for i=1 .. ```).
    please provide the optimized code without including the original C program:
    {search_code}

    """

        else:
            prompt = f"""
    You are a compiler. Previously, we asked you to improve the C program's performance using meaning-preserving transformation:  
    {search_code}
    Here are the different optimized versions you provided, which are ranked(0 being the best): \n"""
            for idx, feedback in enumerate(feedback_list):
                prompt += f"{idx}. Optimized Version:  \n{feedback}\n"
            prompt += """Based on the above rankings, please generate a further optimized version of the code. Put your code in markdown the code block(like ```c for i=1 .. ```).
    Do not define a new function. The original variable does not need to be redefined. If you generate a new variable, please use the double type. 
    """
        messages = [{"role": "user", "content": prompt}]
        response = llm_api[args.llm](messages)
        opt_code = extract_and_concatenate(response)
        print("final_opt_code: \n", opt_code)
        return opt_code, messages


    ########################################################
    def code_after_time_feedback(
        search_code: str, basename: str, feedback_list: list
    ) -> str:
        """
        feedback to llm
        """
        opt_code, messages = get_final_code_from_llm(search_code, feedback_list)
        opt_file = deal_data.code_recombine_lore(
            opt_code,
            basename,
            idx="final",
            org_folder_path=ori_dir,
            opt_folder_path=opt_dir,
        )
        basename_idx = basename + "_final"
        file_path_ = os.path.join(opt_dir, basename_idx + ".after.c")
        error = deal_data.lore_compile(
            file_path=file_path_,
            lore_header_path=ori_dir,
            output_path=os.path.join(opt_dir, basename_idx + ".out"),
        )
        return error, opt_code, messages


    #######################################################
    def first_error_dealing(
        opt_code: str, messages: list, error: str, basename: str, idx: Union[str, int]
    ):
        """
        get new_opt_code from llm,compile it and return error message
        """
        opt_code = get_newoptcode_from_llm(messages, error, opt_code)
        opt_file = deal_data.code_recombine_lore(
            opt_code,
            basename,
            idx=idx,
            org_folder_path=ori_dir,
            opt_folder_path=opt_dir,
        )
        basename_idx = basename + f"_{idx}"
        file_path_ = os.path.join(opt_dir, basename_idx + ".after.c")
        error = deal_data.lore_compile(
            file_path=file_path_,
            lore_header_path=ori_dir,
            output_path=os.path.join(opt_dir, basename_idx + ".out"),
        )
        return error, opt_code


    #######################################################
    def second_error_dealing(
        opt_code: str, messages: list, error: str, basename: str, idx: Union[str, int]
    ):
        """
        get new_opt_code from llm again,compile it and return error message
        """
        opt_code = get_newoptcode_from_llm(messages, error, opt_code, "again")
        opt_file = deal_data.code_recombine_lore(
            opt_code,
            basename,
            idx=idx,
            org_folder_path=ori_dir,
            opt_folder_path=opt_dir,
        )
        basename_idx = basename + f"_{idx}"
        file_path_ = os.path.join(opt_dir, basename_idx + ".after.c")
        error = deal_data.lore_compile(
            file_path=file_path_,
            lore_header_path=ori_dir,
            output_path=os.path.join(opt_dir, basename_idx + ".out"),
        )
        return error, opt_code



    #########################################################
    def data_dealings(
        opt_code: str, idx, error, basename_idx, c1, c2, a, check, run, opts, messages: list
    ):
        if error is None:
            c1[basename_idx] = True
            c2[basename_idx] = True
            exec_path = os.path.join(opt_dir, basename_idx)
            a[basename_idx], check[basename_idx], run[basename_idx] = (
                deal_data.lore_checksum_time_test_exec(
                    exec_path=exec_path + ".out", ori_sum=ori_sum
                )
            )
            if check[basename_idx] == False:
                elemcheck[basename_idx] = False
                run[basename_idx] = float("inf")
            else:
                if deal_data.lore_elemwise(
                    ori_path=os.path.join(ori_dir, basename + ".check.c"),
                    opt_path=os.path.join(opt_dir, basename_idx + ".after.c"),
                ):
                    elemcheck[basename_idx] = True
                else:
                    elemcheck[basename_idx] = False
                    run[basename_idx] = float("inf")

        else:
            error = error[:args.bound]
            c1[basename_idx] = False
            with open("./lore/compile_result/error.txt", "a") as f:
                f.write(
                    "i,basename_idx" + "  " + str(i) + "  " + basename_idx + "\n" + error
                )
            error, opt_code = first_error_dealing(opt_code, messages, error, basename, idx)
            opts[basename_idx] = opt_code
            if error is None:
                c2[basename_idx] = True
                exec_path = os.path.join(opt_dir, basename_idx)
                a[basename_idx], check[basename_idx], run[basename_idx] = (
                    deal_data.lore_checksum_time_test_exec(
                        exec_path=exec_path + ".out", ori_sum=ori_sum
                    )
                )
                if check[basename_idx] == False:
                    elemcheck[basename_idx] = False
                    run[basename_idx] = float("inf")
                else:
                    if deal_data.lore_elemwise(
                        ori_path=os.path.join(ori_dir, basename + ".check.c"),
                        opt_path=os.path.join(opt_dir, basename_idx + ".after.c"),
                    ):
                        elemcheck[basename_idx] = True
                    else:
                        elemcheck[basename_idx] = False
                        run[basename_idx] = float("inf")
            else:
                error = error[:args.bound]
                c2[basename_idx] = False
                a[basename_idx] = False
                check[basename_idx] = False
                elemcheck[basename_idx] = False
                run[basename_idx] = float("inf")
                with open("./lore/compile_result/error.txt", "a") as f:
                    f.write("\n" + error)
                # error, opt_code = second_error_dealing(opt_code,messages,error,basename_idx)
                # if error:
                # error=error[:args.bound]


    #############################################################
    
    with open("./data/raw_data/rag_lore_50_1010.json", "r") as searchcodes:
        searchcodes = json.load(searchcodes)

        for i, data in enumerate(searchcodes):
            search_code, basename = data["code"], data["filename"]
            print(i, basename)
            query_code = data
            ori_sum = data['ori_sum']
            query_code.pop('ori_sum')
            c1, c2, a, check, elemcheck, run, opts = (

                dict(),
                dict(),
                dict(),
                dict(),
                dict(),
                dict(),
                dict(),
            )
            opt_codes=[]
            messages_all=[]
            if args.use_rag:
                results = deal_data.use_rag(query_code,args.es_host,args.index_name)
                # # 打印搜索结果
                # for hit in results["hits"]["hits"]:
                #     print(f"Code: {hit['_source']['code']}")
                #     print(f"Optimized Code: {hit['_source']['opt_code']}")
                #     print(f"Score: {hit['_score']}")
                #     print("-----")
                n = sum(
                        1 for hit in results["hits"]["hits"] if hit["_score"] < args.threshold
                    )
                results_=[results["hits"]["hits"][0],results["hits"]["hits"][1]]
            for j in range(args.method_num):
            ###############################search and prompt making#################################
                if args.use_rag:
                    results_.append(results["hits"]["hits"][j+2])
                    if n == args.search_size:
                        prompt = first_prompt(search_code, use_rag=False)  #########回退机制
                    
                    else:
                        prompt = first_prompt(search_code, args.use_rag, results_)
                    results_.pop()

                else:
                    prompt = first_prompt(search_code, args.use_rag)
                print("prompt:", prompt)
            ###########################get first llm response################################
                messages = [{"role": "user", "content": prompt}]

                response = llm_api[args.llm](messages)

                print("llm_response:\n", response)

                ###################extract opt_code######################################
                opt_code = extract_and_concatenate(response)
                opt_codes.append(opt_code)
                messages_all.append(messages)
            #######################recombine opt_code,return opt_file###############
            for idx, opt_code in enumerate(opt_codes):

                opt_file = deal_data.code_recombine_lore(
                    opt_code,
                    basename,
                    idx,
                    org_folder_path=ori_dir,
                    opt_folder_path=opt_dir,
                )  ########返回opt_dir/ALPBench+ASC+Cortexsuite/1_ALPBench_addMatrixEquals_0.after.c
                ##################compile and return error message#######################
                basename_idx = basename + "_" + str(idx)
                opts[basename_idx] = opt_code
                extension = ".after.c"
                file_path_ = os.path.join(opt_dir, basename_idx + ".after.c")
                error = deal_data.lore_compile(
                    file_path=file_path_,
                    lore_header_path=ori_dir,
                    output_path=os.path.join(opt_dir, basename_idx + ".out"),
                )
                if error:
                    print(error[:args.bound])
                #########################################################################
                try:
                    data_dealings(
                        opt_code,
                        idx,
                        error,
                        basename_idx,
                        c1,
                        c2,
                        a,
                        check,
                        run,
                        opts,
                        messages_all[j],
                    )
                    flag = True
                except Exception as e:
                    logger.exception(f"{basename_idx} datadealing_error:{e}")
                    flag = False
                #########################################################################
            if flag == False:
                continue
            time_feedback = []
            for basename_idx, time_ in run.items():
                if time_ != float("inf"):
                    time_feedback.append((time_, basename_idx))
            time_feedback.sort()
            feedback_list = []
            for time_, basename_idx in time_feedback:
                feedback_list.append(opts[basename_idx])
            error, opt_code, fi_messages = code_after_time_feedback(
                search_code, basename, feedback_list
            )

            opts[basename + "_final"] = opt_code
            try:
                data_dealings(
                    opt_code=opt_code,
                    idx="final",
                    error=error,
                    basename_idx=basename + "_final",
                    c1=c1,
                    c2=c2,
                    a=a,
                    check=check,
                    run=run,
                    opts=opts,
                    messages=fi_messages,
                )
            except Exception as e:
                logger.exception(f"{basename}_final datadealing_error:{e}")
                continue
            temp_rec = {}
            temp_rec["c1"] = c1
            temp_rec["c2"] = c2
            temp_rec["a"] = a
            temp_rec["check"] = check
            temp_rec["elemcheck"] = elemcheck
            temp_rec["run"] = run
            temp_rec["opts"] = opts
            temp_rec["i"] = i

            logger.info(temp_rec)
if __name__ == "__main__":
    main()