generate codes with retrieval or without retrieval

```
python main.py --dataset polybench --method generate_multi --use_rag True 
```

compare results

```
python compare_result.py --dataset polybench 
```

count original checksum and execution time of datasets

```
python count_oritime_and_checksum.py --dataset polybench 
```

count execution time of pluto_optimized versions of datasets

```
python count_plutotime.py  --dataset polybench 
```
Don't forget to fit in your api-key in llm_response.py!

