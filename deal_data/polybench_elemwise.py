import json
import os
import subprocess
import csv

def replace_initial_code(code, new_func, func_name, latter=False):
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

def polybench_exectest_compile(file_path, header_path, polybench_header_path, output_path):
    '''
    !!! make sure your result_path file has been initialized.\n
    example:\n
    polybench_exectest_compile('./polybench/datamining/correlation/correlation.c', './polybench/datamining/correlation', './polybench/utilities', './polybench/datamining/correlation/correlation.out', './polybench/exectest_compile.txt')
    '''
    cmd = ['gcc', '-fopenmp', '-O3', '-fno-tree-dce', '-I', polybench_header_path, '-I', header_path, f'{polybench_header_path}/polybench.c', file_path, '-DPOLYBENCH_DUMP_ARRAYS', '-lm', '-o', output_path]
    print(" ".join(cmd))
    output = subprocess.run(cmd, capture_output=True, text=True)
    if output.returncode:
        return 'error'
    elif 'warning' in output.stderr.lower():
        return 'warning'
    else:
        return 'pass'
        
def polybench2_exectest_exec(exec_path, seed = 1):
    '''
    !!! make sure your result_path file has been initialized.\n
    example:\n
    polybench_exectest_exec('./polybench/datamining/correlation/correlation.out', './polybench/exectest_exec.txt')
    '''
    filename = os.path.basename(exec_path).split(".")[0]
    cmd = [exec_path, f'{seed}']
    # print(" ".join(cmd))
    output = subprocess.run(cmd, capture_output=True, text=True)
    if output.stderr:
        return output.stderr
    
#polybench 单例测试
def check_single_polybench_code(file_path, opt_file_path, include_inf=False):

    filename = os.path.basename(file_path).split('.')[0]
    try:
        polybench_exectest_compile(file_path, os.path.dirname(file_path), 
                                        './polybench/polybench/utilities', f'{os.path.dirname(file_path)}/{filename}.out')
        output = polybench2_exectest_exec(f'{os.path.dirname(file_path)}/{filename}.out')
        
        polybench_exectest_compile(opt_file_path, os.path.dirname(opt_file_path), 
                                    './polybench/polybench/utilities', f'{os.path.dirname(opt_file_path)}/{filename}.out')
        opt_output = polybench2_exectest_exec(f'{os.path.dirname(opt_file_path)}/{filename}.out')
        
        if output == opt_output:
            if not output:
                flag = True
                print(f"{filename} and {filename}_opt are identical.")
            elif output.find("inf") != -1 or output.find("nan") != -1:
                flag = include_inf
                # print(f"{filename} and {filename}_opt have wrong numbers.\n")
            else:
                flag = True
                # print(f"{filename} and {filename}_opt are identical.")
        else:
            flag = False
            # print(f"{filename} and {filename}_opt are different.")
        return flag
        
    except FileNotFoundError as e:
        print("语法错误")
        return False
    
    except Exception as e:
        print(f"未知错误: {e}")
        return False
    
def polybench_elemwise(ori_path, opt_path):
    # ori_path 原来数据集的path ./s000.check.c
    # opt_path 新的path        ./s000_0.after.c
    basename = os.path.basename(ori_path)
    filename = basename.split('.')[0]
    ori_dirname = os.path.dirname(ori_path)
    opt_dirname = os.path.dirname(opt_path)
    input_path = './polybench_init_25.json'
    with open(input_path, 'r') as f1, open(ori_path, 'r') as f2, open(opt_path, 'r') as f3:
        json_data = json.loads(f1.read())
        ori_code = f2.read()
        opt_code = f3.read()

    code_list = json_data[basename]['init_func']
    # k = len(code_list)
    k = 5
    for code in code_list[:k]:
        init_code, _, _ = extract_initial_code(code, "{", "}", "static\nvoid init_array")
        # 替换代码
        if init_code is not None:
            new_ori_code = replace_initial_code(ori_code, code, "static\nvoid init_array") 
            new_opt_code = replace_initial_code(opt_code, code, "static\nvoid init_array")
        else:
            init_code, _, _ = extract_initial_code(code, "{", "}", "static void init_array")
            new_ori_code = replace_initial_code(ori_code, code, "static void init_array") 
            new_opt_code = replace_initial_code(opt_code, code, "static\nvoid init_array")
        
        new_opt_path = os.path.join(opt_dirname, "polybench_opt_temp.c")
        new_ori_path = os.path.join(ori_dirname, "polybench_temp.c")

        with open(new_ori_path, 'w') as f1, open(new_opt_path, 'w') as f2:
            f1.write(new_ori_code)
            f2.write(new_opt_code)

        flag = check_single_polybench_code(new_ori_path, new_opt_path, include_inf=True)
        # os.remove(new_ori_path)
        # os.remove(new_opt_path)
        if not flag:
            return False 
    return True

if __name__ == '__main__':
    flag = polybench_elemwise("./polybench/polybench/linear-algebra/blas/syr2k/syr2k.check.c", "./polybench/polybench/linear-algebra/blas/syr2k/syr2k_final.after.c")
    print(flag)
