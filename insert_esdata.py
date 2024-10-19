import json,argparse
from elasticsearch import Elasticsearch, helpers  

def delete_index_if_exists(es, index_name):  
    if es.indices.exists(index=index_name):  
        es.indices.delete(index=index_name)  
        print(f"Deleted index '{index_name}'")  
    else:  
        print(f"Index '{index_name}' does not exist.")  

def create_index(es, index_name, mapping):  
    es.indices.create(index=index_name, body=mapping)  
    print(f"Created index '{index_name}'")  

def read_and_store_json(file_path1, es, index_name):  
    with open(file_path1, "r") as f:  
        datas = json.load(f)  
        actions = [  
            {  
                "_index": index_name,  
                "_id": doc.get("filename"),  # 使用 filename 作为文档ID  
                "_source": doc,  
            }  
            for doc in datas[:50000]  
        ]  
        actions1 = [  
            {  
                "_index": index_name,  
                "_id": doc.get("filename"),  # 使用 filename 作为文档ID  
                "_source": doc,  
            }  
            for doc in datas[50000:]  
        ]  
    helpers.bulk(es, actions)
    helpers.bulk(es, actions1)    
    print(f"Inserted {len(datas)} documents into the index '{index_name}'.")  

def main(args):  
    es_host = args.es_host  # 替换为您的Elasticsearch服务器地址  
    index_name = args.index_name 
    mapping_file = args.mapping_file
    file_path1 = args.file_path

    es = Elasticsearch([es_host])  

    # 1. 删除现有索引  
    delete_index_if_exists(es, index_name)  

    # 2. 创建新的索引并应用mapping  
    with open(mapping_file, "r") as f:  
        mapping = json.load(f)  
        create_index(es, index_name, mapping)  

    # 3. 插入数据  
    read_and_store_json(file_path1, es, index_name)  

if __name__ == "__main__":  
    argparser= argparse.ArgumentParser()
    argparser.add_argument("--file_path",type=str,default="./data/raw_data/rag_189487_1014.json")
    argparser.add_argument("--es_host",type=str,default="http://10.214.241.234:9210")
    argparser.add_argument("--index_name",type=str,default="rag_189487_1014")
    argparser.add_argument("--mapping_file",type=str,default="./settings/mappings/mapping7.json")
    args=argparser.parse_args()
    
    main(args)