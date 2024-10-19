import re
############################################
def extract_code_and_text(text):
# 匹配文本说明
    text_pattern = re.compile(r'###([\s\S]*?)```c')
    text_blocks = text_pattern.findall(text)

    # 匹配代码块
    code_pattern = re.compile(r'```c([\s\S]*?)```')
    code_blocks = code_pattern.findall(text)


    print("Text Blocks:")
    print(text_blocks)

    print("\nCode Blocks:")
    print(code_blocks)
    codes=[]
    for code,text in zip(code_blocks,text_blocks):
        codes.append("/*"+text+"*/"+code)
    return codes

def extract_and_concatenate(text):
    # 正则表达式匹配 ```c 和 ``` 包围的内容
    pattern = r'```c(.*?)```'
    
    # 提取匹配的内容
    code_blocks = re.findall(pattern, text, re.DOTALL)
    
    
    # 使用正则表达式替换提取的部分为空字符串
    remaining_text = re.sub(pattern, '', text, flags=re.DOTALL)
    
    # 拼接剩余的内容并移除多余空白
    remaining_text = remaining_text.strip()

    
    code= "/*"+remaining_text+"*/\n"+code_blocks[0]
    print("final_response:\n",code)
    return code
