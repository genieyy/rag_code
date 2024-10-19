import nltk, json
from nltk.translate.bleu_score import sentence_bleu, SmoothingFunction
from zhengze import  zhengze

# 确保已经下载了nltk的数据
# nltk.download('punkt')
def bleu(candidate, reference):

    # 参考翻译
    reference = [
        reference.split(),
    ]

    # 生成翻译
    candidate = candidate.split()

    # 计算BLEU分数
    smooth = SmoothingFunction().method4
    score = sentence_bleu(reference, candidate, smoothing_function=smooth)

    print(f"BLEU score: {score}")
    return score

if __name__ == '__main__':
    with open("./sampled_output.json") as f:
        data = json.load(f)
        dataall = {}

        references = []
        candidates = []
        rag_code = []
        norag_code = []
        for i in data:
            references.append(i["code"])
            candidates.append(i["opt_code"])
            basename = i["filename"]
            with open(
                f"./loop_transformation_classifier_pluto_code/{basename}.after.c",
                "r",
            ) as f1, open(
                f"./loop_transformation_classifier_pluto_code_no_rag/{basename}.after.c",
                "r",
            ) as f2:
                rag_code.append(zhengze(f1.read()))
                norag_code.append(zhengze(f2.read()))

        dataall["ragcode"] = rag_code
        dataall["noragcode"] = norag_code
        dataall["references"] = references
        dataall["candidates"] = candidates
        n=len(dataall["references"])
        bleurag_re = 0
        bleunorag_re = 0
        bleurag_ca= 0
        bleunorag_ca= 0
        for i in range(n):
            bleurag_re+=bleu(dataall["references"][i],dataall["ragcode"][i])
            bleunorag_re+=bleu(dataall["references"][i],dataall["noragcode"][i])
            bleurag_ca+=bleu(dataall["candidates"][i],dataall["ragcode"][i])
            bleunorag_ca+=bleu(dataall["candidates"][i],dataall["noragcode"][i])
        print(bleurag_re/n,bleunorag_re/n,bleurag_ca/n,bleunorag_ca/n)
        
        
