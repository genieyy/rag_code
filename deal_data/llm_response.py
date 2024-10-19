from openai import OpenAI
import os

from transformers import AutoTokenizer, AutoModelForCausalLM


def get_response_deepseek(new_massages):
    client = OpenAI(api_key="sk-0c621ec78f624e818dcd110f5a12993d", base_url="https://api.deepseek.com/beta")

    response = client.chat.completions.create(
        model="deepseek-coder",
        messages=[
            {"role": "system", "content": "You are a helpful compiler"},
        ]+new_massages,
        stream=False,
          
    )
    return response.choices[0].message.content

def get_response_gpt4o(new_massages):
    client = OpenAI(api_key="sk-RNKvKjo6bwNpGC79E3D904Df5f32464fA08b251235782aFd",base_url="http://29qg.com/v1/")
    response = client.chat.completions.create(
        model="gpt-4o-2024-08-06",
        messages=[
            {"role": "system", "content": "You are a helpful compiler"},
        ]+new_massages,
        stream=False,
          
    )
    return response.choices[0].message.content





llm_api={
    "deepseek-coder":get_response_deepseek,
    "gpt-4o":get_response_gpt4o

}

