import subprocess
import os
import csv
import re
import json, math

def tsvc_compile(file_path, output_path):
    """
    example:\n
    tsvc__compile('./new_filter/s000.c', './new_filter/s000.out')
    """
    cmd = [
        "gcc",
        "-fopenmp",
        "-O3",
        "-I",
        os.path.join(tsvc_dir, "cfiles"),
        os.path.join(tsvc_dir, "cfiles/dummy.c"),
        os.path.join(tsvc_dir, "cfiles/common.c"),
        file_path,
        "-lm",
        "-o",
        output_path,
    ]
    output = subprocess.run(cmd, capture_output=True, text=True)

    if output.returncode:
        return output.stderr
    else:
        return None


def tsvc_exectest_exec(exec_path):
    """
    !!! make sure your result_path file has been initialized.\n
    example:\n
    tsvc_exectest_exec('./new_filter/s000.out')
    """

    import subprocess

    cmd = f"{exec_path}"

    try:
        result = subprocess.run(cmd, capture_output=True, timeout=200)

        if result.returncode == 0:
            return True
        else:
            return False
    except subprocess.TimeoutExpired:
        return "timeout"


def tsvc_checksum_time_test_exec(exec_path, ori_sum=None, plutotime=None):
    """
    !!! make sure your result_path file has been initialized.\n
    example:\n
    tsvc_checksumtest_exec('./new_filter/s000.checksum')
    """

    cmd = [exec_path]
    output = subprocess.run(cmd, capture_output=True, text=True)
    time = float(output.stdout.split("\n")[1].split("\t")[1])
    if plutotime:
        time = float(output.stdout.split("\n")[1].split("\t")[1])
        return time
    if ori_sum is None:  ###############没有ori_sum,则用于计算ori_sum
        return float(output.stdout.split("\n")[1].split("\t")[-1]), time

    if not output.returncode:
        sum = float(output.stdout.split("\n")[1].split("\t")[-1])
 
        if sum != ori_sum:
            if not math.isclose(sum, ori_sum):
                return False, float("inf")
        return True, time
    else:
        return False, float("inf")

def polybench_compile(file_path, header_path, polybench_header_path, output_path):
    """
    !!! make sure your error_path file has been initialized.\n
    example:\n
    polybench_checksumtest_compile(
    './polybench/datamining/correlation/correlation.c',
    './polybench/datamining/correlation',
    './polybench/utilities',
    './polybench/datamining/correlation/correlation.checksum')
    """
    file=os.path.basename(file_path).split(".")[0].split("_")[0]
    
    extra_large = [
        "atax",
        "bicg",
        "mvt",
        "gemm",
        "gemver",
        "gesummv",
        "durbin",
        "trisolv",
        "jacobi-1d",
    ]
    if file in extra_large:
        cmd = [
            "gcc",
            "-fopenmp",
            "-O3",
            "-I",
            polybench_header_path,
            "-I",
            header_path,
            f"{polybench_header_path}/polybench.c",
            file_path,
            # "-DPOLYBENCH_DUMP_ARRAYS",
            "-DPOLYBENCH_TIME",
            "-DPOLYBENCH_CHECKSUM_ARRAYS",
            "-DEXTRALARGE_DATASET"
            "-lm",
            "-o",
            output_path,
        ]
    else:
        cmd = [
            "gcc",
            "-fopenmp",
            "-O3",
            "-I",
            polybench_header_path,
            "-I",
            header_path,
            f"{polybench_header_path}/polybench.c",
            file_path,
            # "-DPOLYBENCH_DUMP_ARRAYS",
            "-DPOLYBENCH_TIME",
            "-DPOLYBENCH_CHECKSUM_ARRAYS",
            "-lm",
            "-o",
            output_path,
        ]
    output = subprocess.run(cmd, capture_output=True, text=True)
    if output.returncode:
        return output.stderr
    elif "warning" in output.stderr.lower():
        return None
    else:
        return None


def polybench_checksum_time_test_exec(
    exec_path, ori_sum=None, plutotime=False, seed=1
):  # -> list | tuple[Literal[False], float] | tuple[Literal[True]...:# -> list | tuple[Literal[False], float] | tuple[Literal[True]...:# -> list | tuple[Literal[False], float] | tuple[Literal[True]...:# -> list | tuple[Literal[False], float] | tuple[Literal[True]...:# -> list | tuple[Literal[False], float] | tuple[Literal[True]...:
    """
    包括能否执行，能否通过checksum,执行时间
    !!! make sure your result_path file has been initialized.\n
    example:\n
    polybench_checksumtest_exec('./polybench/datamining/correlation/correlation.checksum')
    """
    cmd = [exec_path]
    
    try:
        output = subprocess.run(cmd, capture_output=True, text=True, timeout=100)

        if output.returncode:
            return False, False, float("inf")
    except subprocess.TimeoutExpired:
        return float("inf")

   
    sum = re.findall(r"-?\d+\.\d+|-?\d+", output.stderr)
    time = float(output.stdout.split("\n")[0])
    result = []
    for s in sum:
        result.append(float(s))
    if plutotime:
        return time
    elif ori_sum is None:
        return result,time

    else:
        for i in range(len(result)):
            if result[i] != ori_sum[i]:
                if not math.isclose(result[i], ori_sum[i]):
                    return True, False, float("inf")
        return True, True, time


def lore_compile(file_path, lore_header_path, output_path):
    """
    能否正确编译，不能返回error
    """
    cmd = [
        "gcc",
        "-O3",
        "-fopenmp",
        "-I",
        lore_header_path,
        os.path.join(lore_header_path, "lore.c"),
        file_path,
        "-lm",
        "-DLORE_CHECKSUM_ARRAYS",
        #  "-DLORE_DUMP_ARRAYS",
        "-o",
        output_path,
    ]
    output = subprocess.run(cmd, capture_output=True, text=True)

    if output.returncode:
        return output.stderr
    elif "warning" in output.stderr.lower():
        return None
    else:
        return None


def lore_checksum_time_test_exec(
    exec_path, ori_sum=None, plutotime=False, seed=1
):  # -> list | tuple[Literal[False], float] | tuple[Literal[True]...:# -> list | tuple[Literal[False], float] | tuple[Literal[True]...:# -> list | tuple[Literal[False], float] | tuple[Literal[True]...:# -> list | tuple[Literal[False], float] | tuple[Literal[True]...:# -> list | tuple[Literal[]]]]]]:
    """
    包括能否执行，能否通过checksum,执行时间
    """
    cmd = [exec_path]
    try:
        output = subprocess.run(cmd, capture_output=True, text=True,timeout=100)

        if output.returncode:
            return False, False, float("inf")
    except subprocess.TimeoutExpired:
        return float("inf")
   
    sum = re.findall(r"-?\d+\.\d+|-?\d+", output.stderr)
    time = float(output.stdout.split("\n")[0])
    result = []
    for s in sum:
        result.append(float(s))
    if plutotime:
        return time
    elif ori_sum is None:
        return result,time

    else:
        for i in range(len(result)):
            if result[i] != ori_sum[i]:
                if not math.isclose(result[i], ori_sum[i]):
                    return True, False, float("inf")
        return True, True, time



if __name__ == "__main__":

    tsvc_dir="./tsvc"
    polybench_dir="./polybench"
    lore_compile(
        "lore/LORE_artificial/Freebench/1_Freebench_pifft4.check.c",
        "lore/LORE_artificial",
        "./lore.out",
    )
    lore_checksum_time_test_exec("./lore.out")
