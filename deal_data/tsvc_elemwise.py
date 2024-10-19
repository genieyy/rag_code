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
                    if file.endswith(".c") and not file.endswith("_code.c") and not file.endswith("tsvc_temp.c"):
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

def tsvc2_exectest_compile(file_path, header_path, output_path, common_path="common.c", args=None):
    cmd = ['gcc', '-fopenmp', '-O3', '-I', header_path, '-I', header_path, f'{common_path}', 
           file_path, f'{header_path}/dummy.c', "-DTSVC_DUMP_ARRAYS",'-lm', '-o', output_path]
    if args:
        cmd.extend(args)
    print(" ".join(cmd))
    output = subprocess.run(cmd, capture_output=True, text=True)
    if output.returncode:
        return 'error'
    elif 'warning' in output.stderr.lower():
        return 'warning'
    else:
        return 'pass'
        
def tsvc2_exectest_exec(exec_path):
    cmd = [f"{exec_path}"]
    print(" ".join(cmd))
    output = subprocess.run(cmd, capture_output=True, text=True)
    if output.stderr:
        return output.stderr

def compare_by_relative_tolerance(output, output_opt):
    if not(output and output_opt):
        return False######如果其中一个为None，返回False
    output_list = output.strip().split(" ")
    opt_output_list = output_opt.strip().split(" ")
    for num1, num2 in zip(output_list, opt_output_list):
        if num1 == "nan" or num1 == "inf": continue
        flag = math.isclose(float(num1), float(num2))
        if not flag:
            return False
    return True

# tsvc 单例测试
def check_single_tsvc_code(filename="s112", common_path="common.c", ori_path="", opt_path=""):
    try: 
        tsvc2_exectest_compile(ori_path, './tsvc/cfiles',                                    
            './tsvc_temp.out', common_path=common_path)
        output = tsvc2_exectest_exec('./tsvc_temp.out')
        os.remove('./tsvc_temp.out')
        tsvc2_exectest_compile(opt_path, './tsvc/cfiles',                                 
            './tsvc_temp_opt.out', common_path=common_path)
        output_opt = tsvc2_exectest_exec('./tsvc_temp_opt.out')
        os.remove('./tsvc_temp_opt.out')

        if output == output_opt:
            if not output:
                flag = True
                print(f"{filename} and {filename}_opt are identical.")
            elif output.find("inf") != -1 or output.find("nan") != -1:
                flag = True
                print(f"{filename} and {filename}_opt have wrong numbers.\n")
            else:
                flag = True
                print(f"{filename} and {filename}_opt are identical.")
        else:
            flag = compare_by_relative_tolerance(output, output_opt)
            if flag:
                print(f"{filename} and {filename}_opt are identical.")
            else:
                print(f"{filename} and {filename}_opt are different.")
        return flag
    except FileNotFoundError as e:
        print(f"语法错误: {e}")
        return False
    except Exception as e:
        print(f"未知错误: {e}")
        return False
    
def tsvc_elemwise(ori_path, opt_path):
    # ori_path 原来数据集的path ./s000.check.c
    # opt_path 新的path        ./s000_0.after.c
    basename = os.path.basename(ori_path)
    filename = basename.split('.')[0]
    input_path = '.data/elemwise/set1d_init.json'
    template_path = './tsvc/cfiles/common_template.c'
    num_pass = 0
    with open(input_path, 'r') as f1, open(template_path, 'r') as f2:
        json_data = json.loads(f1.read())
        template = f2.read()
    code_list = json_data[basename]
    # k = len(code_list)
    k = 5
    for code in code_list[:k]:
        new_code = replace_tsvc_code(template, code, "void set_1d_array", latter=True) 
        temp_path = "./tsvc_temp.c"
        with open(temp_path, 'w') as f:
            f.write(new_code)
        flag = check_single_tsvc_code(filename, common_path=temp_path, ori_path=ori_path, opt_path=opt_path)
        # 如果测试通过就删除，否则保留以供查看
        if not flag:
            return False 
        os.remove(temp_path)
    return True

if __name__ == '__main__':
    ori_path = "./tsvc/cfiles/s000.check.c"
    opt_path = "./tsvc/opt_files/s000_0.after.c"
    flag = tsvc_elemwise(ori_path, opt_path)
    print(flag)
    