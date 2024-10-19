import json
import os
import subprocess
import csv
import math


def list_files_in_directory(directory_path, include_dir=False, suffix=None, exclude_suffix=None, traverse_deeper=True, words=None):
    c_files = []
    try:
        for root, dirs, files in os.walk(directory_path):
            for file in files:
                if suffix is not None and not file.endswith(suffix): continue
                if exclude_suffix is not None and file.endswith(exclude_suffix): continue
                if words:
                    if file.find(words) != -1:
                        if include_dir:
                            c_files.append(os.path.join(root, file))
                        else:
                            c_files.append(file)
                else:
                    if file.endswith(".c") and not file.endswith("_code.c") and not file.endswith("temp.c"):
                        if include_dir:
                            c_files.append(os.path.join(root, file))
                        else:
                            c_files.append(file)
            # 根据布尔变量决定是否遍历更深的层
            if not traverse_deeper:
                # 只遍历当前目录，因此打破循环，不进入子目录
                break
        # 按文件名排序
        c_files.sort()
    except FileNotFoundError as e:
        print(f"Error: {e}")
    return c_files

def replace_tsvc_code(code, new_func, func_name, latter=False):
    # 寻找tsvc_init_array函数
    _, func_start, func_end = extract_initial_code(code, "{", "}", func_name, include_symbol=True, latter=latter)
    return code[0:func_start] + new_func + code[func_end+1:]

def extract_initial_code(code, start_symbol, end_symbol, entry="static void init_array", 
                         include_symbol=False, latter=False):
    start_index = code.find(entry)
    if latter:
        start_index = code.find(entry, start_index + 1)
    num = 0
    flag = False
    for i in range(start_index, len(code)):
        if code[i] == start_symbol:
            num += 1
            if not flag: 
                start_symbol_index = i
            flag = True
        elif code[i] == end_symbol:
            num -= 1
        if flag and num == 0:
            end_index = i
            end_symbol_index = end_index
            if include_symbol:
                return code[start_index: end_index+1], start_index, end_index
            else:
                return code[start_index: end_index+1], start_symbol_index, end_symbol_index
    return None, None, None

def lore_exectest_compile(file_path, lore_header_path, output_path, common_path):

    cmd = ['gcc', '-fopenmp',  '-O3', '-I', lore_header_path, f'{common_path}', file_path, '-lm', '-DLORE_DUMP_ARRAYS', '-o', output_path]
    print(" ".join(cmd))
    output = subprocess.run(cmd, capture_output=True, text=True)
    if output.returncode:
        return 'error'
    elif 'warning' in output.stderr.lower():
        return 'warning'
    else:
        return 'pass'
        
def lore_exectest_exec(exec_path):
    cmd = [exec_path]
    output = subprocess.run(cmd, capture_output=True, text=True)
    if output.stderr:
        return output.stderr

def compare_by_relative_tolerance(output, output_opt):
    if not(output and output_opt):
        return False######如果其中一个为None，返回False
    output_list = output.strip().split(" ")
    opt_output_list = output_opt.strip().split(" ")
    import time
    start = time.time()
    for num1, num2 in zip(output_list, opt_output_list):
        flag = math.isclose(float(num1), float(num2))
        if not flag:
            return False
    return True
    
# tsvc 单例测试
def check_single_lore_code(ori_path="", opt_path="", common_path="common.c"):
    try: 
        lore_exectest_compile(ori_path, './lore/LORE_artificial', './temp.out', common_path=common_path)
        output = lore_exectest_exec('./temp.out')
        os.remove('./temp.out')
        lore_exectest_compile(opt_path, './lore/LORE_artificial', './temp_opt.out', common_path=common_path)
        output_opt = lore_exectest_exec('./temp_opt.out')
        os.remove('./temp_opt.out')
        if output == output_opt:
            if not output:
                flag = True
                print(f"identical.")
            if output.find("inf") != -1 or output.find("nan") != -1:
                flag = True
                print(f"Wrong numbers.\n")
            else:
                flag = True
                print(f"Identical.")
        else:
            flag = compare_by_relative_tolerance(output, output_opt)
            if flag: print(f"Identical.")
            else: print(f"Different.")
        return flag
    except FileNotFoundError as e:
        print(f"语法错误: {e}")
        return False
    except Exception as e:
        print(f"未知错误: {e}")
        return False
    
def lore_elemwise(ori_path, opt_path):
    # ori_path 原来数据集的path ./s000.check.c
    # opt_path 新的path        ./s000_0.after.c
    basename = os.path.basename(ori_path)
    filename = basename.split('.')[0]
    input_path = './data/elemwise/lore_init_25.json'
    template_path = './data/elemwise/lore_template.c'
    with open(input_path, 'r') as f1, open(template_path, 'r') as f2:
        json_data = json.loads(f1.read())
        template = f2.read()
    code_list = json_data[basename]['code']
    # k = len(code_list)
    k = 10
    for code in code_list[:k]:
        new_code = replace_tsvc_code(template, code, "void init_array_1d") 
        temp_path = "./temp.c"
        with open(temp_path, 'w') as f:
            f.write(new_code)
        flag = check_single_lore_code(ori_path=ori_path, opt_path=opt_path, common_path=temp_path)
        # 如果测试通过就删除，否则保留以供查看
        if not flag:
            return False 
        os.remove(temp_path)
    return True

if __name__ == '__main__':
    # file_list = list_files_in_directory("code", include_dir=True, suffix=".check.c")
    # with open("./output/pluto-2024-0816.txt", "w") as f: pass
    # for ori_path in file_list:
    #     print(ori_path)
    #     opt_path = os.path.join("pluto_code", os.path.basename(ori_path).replace(".check.c", ".pluto.c"))
    #     if not os.path.exists(opt_path):
    #         continue
    #     flag = elemwise(ori_path, opt_path)
    #     with open("./output/pluto-2024-0816.txt", "a") as f:
    #         if flag:
    #             f.write(f"{ori_path} and {opt_path} are identical.\n")
    #         else:
    #             f.write(f"{ori_path} and {opt_path} are different.\n")
    #     print(flag)
    
    ori_path = "./lore/LORE_artificial/ALPBench+ASC+Cortexsuite/1_ALPBench_subtractClassMean.check.c"
    opt_path = "./lore/LORE_artificial/ALPBench+ASC+Cortexsuite/1_ALPBench_subtractClassMean_4.after.c"
    flag = lore_elemwise(ori_path, opt_path)
    print(flag)
    