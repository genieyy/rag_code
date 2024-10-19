import insert_esdata, json
from elasticsearch import Elasticsearch, helpers


def use_rag(query_code,es_host="http://10.214.241.234:9210",index_name="rag_99125_1010") -> bool:


    es = Elasticsearch([es_host])

    
    search_query = {
        "size": 5,  # 设置返回的结果数
        # "explain": True,
        "query": {
            "function_score": {
                "query": {
                    "bool": {
                        "must": [
                            {
                                "match": {
                                    "code": {
                                        "query": query_code["code"],
                                        "boost": 1,  # BM25基本分数
                                    }
                                }
                            }
                        ]
                    }
                },
                "functions": [
                    {
                        "script_score": {
                            "script": {
                                "source": """  
                                double score = 0;  // 使用BM25的基本分数  
                
                                // 获取查询的info字段的长度
                                int queryInfoSize = params['info'].size();

            

                                // 遍历查询的info字段
                                for (int i = 0; i < queryInfoSize; i++) {
                                    Map queryInfoItem = params['info'][i];

                                    List queryScheduleConst = queryInfoItem['schedule_const'];
                                    List queryScheduleVar = queryInfoItem['schedule_iterator'];
                                    List queryIndexesConstRead = queryInfoItem['indexes_const_read'];
                                    List queryIndexesNonzeroVarRead= queryInfoItem['indexes_nonzero_iterator_read'];
                                    List queryIndexesNonzeroVarWrite= queryInfoItem['indexes_nonzero_iterator_write'];
                                    List queryIndexesConstWrite= queryInfoItem['indexes_const_write'];
                                    
                                    // 获取文档的info字段的长度
                                    int docInfoSize = params._source['info'].size();

                                    // 遍历文档中的info字段
                                    if(i<docInfoSize) {
                                        Map docInfoItem = params._source['info'][i];

                                        List docScheduleConst = docInfoItem['schedule_const'];
                                        List docScheduleVar = docInfoItem['schedule_iterator'];
                                        List docIndexesConstRead = docInfoItem['indexes_const_read'];
                                        List docIndexesNonzeroVarWrite = docInfoItem['indexes_nonzero_iterator_write'];
                                        List docIndexesNonzeroVarRead = docInfoItem['indexes_nonzero_iterator_read'];
                                        List docIndexesConstWrite = docInfoItem['indexes_const_write'];
                                        

                                        // schedule_const重叠加分
                                        Set overlapSchedulesConst = new HashSet(queryScheduleConst);
                                        overlapSchedulesConst.retainAll(docScheduleConst);
                                        if(!overlapSchedulesConst.isEmpty()){
                                            
                                            score += 2 * overlapSchedulesConst.size()/docScheduleConst.size();
                                        }

                                        // schedule_iterator重叠加分
                                        Set overlapSchedulesVar = new HashSet(queryScheduleVar);
                                        overlapSchedulesVar.retainAll(docScheduleVar);
                                        if (!overlapSchedulesVar.isEmpty()) {
                                            
                                            score += 3 * overlapSchedulesVar.size()/docScheduleVar.size();
                                        }
                                        
                                        // indexes_const_read重叠加分
                                        Set overlapIndexesConstRead = new HashSet(queryIndexesConstRead);
                                        overlapIndexesConstRead.retainAll(docIndexesConstRead);
                                        if (!overlapIndexesConstRead.isEmpty()){
                                            score += 4 * overlapIndexesConstRead.size()/docIndexesConstRead.size();
                                            }
                                        
                                        // indexes_nonzero_iterator_read重叠加分
                                        Set overlapIndexesNonzeroVarRead = new HashSet(queryIndexesNonzeroVarRead);
                                        overlapIndexesNonzeroVarRead.retainAll(docIndexesNonzeroVarRead);
                                        if (!overlapIndexesNonzeroVarRead.isEmpty()){
                                            score += 4 * overlapIndexesNonzeroVarRead.size()/docIndexesNonzeroVarRead.size();
                                        }
                                        
                                        
                                        // indexes_const_write重叠加分
                                        Set overlapIndexesConstWrite = new HashSet(queryIndexesConstWrite);
                                        overlapIndexesConstWrite.retainAll(docIndexesConstWrite);
                                        if (! overlapIndexesConstWrite.isEmpty()){
                                            score += 6 * overlapIndexesConstWrite.size()/docIndexesConstWrite.size();
                                        
                                        }
                                        
                                        // indexes_nonzero_iterator_write重叠加分
                                        Set overlapIndexesNonzeroVarWrite = new HashSet(queryIndexesNonzeroVarWrite);
                                        overlapIndexesNonzeroVarWrite.retainAll(docIndexesNonzeroVarWrite);
                                        if (!overlapIndexesNonzeroVarWrite.isEmpty()){
                                            score += 6 * overlapIndexesNonzeroVarWrite.size()/docIndexesNonzeroVarWrite.size();
                                        
                                        }
                                        
                                    }    
                                } 

                                return score*0.6+_score;

                            """,
                                "params": query_code,
                            }
                        }
                    }
                ],
                "boost_mode": "replace",
            }
        },
    }

    # print(json.dumps(search_query, indent=2, ensure_ascii=False))


    res = es.search(index=index_name, body=search_query,request_timeout=30)

    return res
        
