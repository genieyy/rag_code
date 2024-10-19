import argparse

if __name__ == "__main__":
    argparser = argparse.ArgumentParser()
    argparser.add_argument(
        "--dataset",
        type=str,
        default="polybench",
        choices=["polybench", "lore", "tsvc"],
    )
    argparser.add_argument(
        "--method",
        type=str,
        default="multi_generate",
        choices=["generate_multi", "multi_generate"],
    )
    argparser.add_argument(
        "--search_size", type=int, default=5, help="number of retrieval results"
    )
    argparser.add_argument(
        "--threshold", type=float, default=15, help="threshold for RAG"
    )
    argparser.add_argument(
        "--method_num", type=int, default=3, help="number of methods to generate"
    )
    argparser.add_argument(
        "--bound",
        type=int,
        default=4096,
        help="max length of compiling error returned to llm",
    )
    argparser.add_argument("--use_rag", type=bool, default=True, help="use rag or not")
    argparser.add_argument(
        "--llm",
        type=str,
        default="deepseek-coder",
        help="choose llm",
        choices=["deepseek-coder", "gpt-4o"],
    )
    argparser.add_argument(
        "--es_host", type=str, default="http://10.214.241.234:9210"
    )
    argparser.add_argument("-i", "--index_name", type=str,default="rag_189487_1014")
    args = argparser.parse_args()
    if args.dataset == "polybench":
        if args.method == "generate_multi":
            from generate_multi_polybench import main
        elif args.method == "multi_generate":
            from multi_generate_polybench import main
    elif args.dataset == "lore":
        if args.method == "generate_multi":
            from generate_multi_lore import main
        elif args.method == "multi_generate":
            from multi_generate_lore import main
    else:
        if args.method == "generate_multi":
            from generate_multi_tsvc import main
        else:
            from multi_generate_tsvc import main
    main(args)
