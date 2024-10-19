import subprocess
import json

with open("benchmark_list", 'r') as f:
    files = f.read().split()

contents = []
errors = []

for file in files:
    output_path = f'./{file[:-2]}.check.out'
    compile_cmd = ['gcc', '-O3', '-fopenmp', '-ffast-math', '-I', '.', './lore.c', f'./{file[:-2]}.check.c', '-lm', '-DLORE_CHECKSUM_ARRAYS', '-o', output_path]
    # compile_cmd = ['clang', '-O3', '-fopenmp', '-ffast-math', '-I', '.', './lore.c', f'./{file[:-2]}.check.c', '-lm', '-o', output_path]
    
    compile_output = subprocess.run(compile_cmd, capture_output=True, text=True)
    
    if compile_output.returncode:
        errors.append(f'compile error:{file}\n')
    else:
        exec_output = subprocess.run([output_path], capture_output=True, text=True)
        
        if exec_output.returncode:
            errors.append(f'exec error:{file}\n')
        else:
            contents.append({'file': file, 'time': exec_output.stdout, 'sum': exec_output.stderr}) # 'sum': exec_output.stderr

with open("error.txt", 'w') as f:
    f.writelines(errors)
         
with open("output.jsonl", 'w') as f:
    for l in contents:
        f.write(json.dumps(l) + '\n')
