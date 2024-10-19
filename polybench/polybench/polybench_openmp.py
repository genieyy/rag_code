import os
import re
import subprocess
import json
from dotenv import load_dotenv

compiler = 'gcc'

if compiler == 'gcc' or 'clang':
    compiler_cmd = [compiler, '-fopenmp', '-O3']
elif compiler == 'icx':
    load_dotenv()
    os.environ['PATH'] = '/opt/intel/oneapi/compiler/2022.1.0/linux/bin:' + os.environ['PATH']
    compiler_cmd = [compiler, '-qopenmp', '-O3', '-xHost']


with open("utilities/benchmark_list", 'r') as f:
    files = f.read().split()

contents = []
errors = []

for file in files:
    header_path = '/'.join(file.split('/')[1:-1])
    
    output_path = f'{file[:-2]}.out'
    
    compile_cmd = compiler_cmd + ['-I', 'utilities', '-I', header_path, 'utilities/polybench.c', f'{file[:-2]}.check.c', '-DPOLYBENCH_TIME', '-DPOLYBENCH_CHECKSUM_ARRAYS', '-lm', '-o', output_path]
    
    compile_output = subprocess.run(compile_cmd, capture_output=True, text=True)
    
    if compile_output.returncode:
        # print('compile: ', file)
        # print(" ".join(compile_cmd))
        errors.append(f'compile error:{file}\n')
    else:
        exec_cmd = [output_path]
        
        exec_output = subprocess.run(exec_cmd, capture_output=True, text=True)
        
        if exec_output.returncode:
            # print('exec: ', file)
            errors.append(f'exec error:{file}\n')
        else:
            sum = re.findall(r'begin dump:(((?!begin dump).)*)\nend   dump', exec_output.stderr, re.DOTALL)
            contents.append({'code': file, 'time': exec_output.stdout[:-1], 'checksum': [s[0].split() for s in sum]})
            
            # for _ in range(2):
            #     exec_output = subprocess.run(exec_cmd, capture_output=True, text=True)
            #     sum = re.findall(r'begin dump:(((?!begin dump).)*)\nend   dump', exec_output.stderr, re.DOTALL)
            #     contents.append({'code': file, 'time': exec_output.stdout[:-1], 'checksum': [s[0].split() for s in sum]})

with open("error.txt", 'w') as f:
    f.writelines(errors)
         
with open("output.jsonl", 'w') as f:
    for l in contents:
        f.write(json.dumps(l) + '\n')