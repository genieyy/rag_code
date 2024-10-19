import re
def zhengze(code):

    pattern = re.compile(r'#pragma\s+scop\s*(.*?)\s*#pragma\s+endscop', re.DOTALL)

# 查找匹配的内容
    matches = pattern.findall(code)

# 输出匹配结果
    return matches[0].strip()