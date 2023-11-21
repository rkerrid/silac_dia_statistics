# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 12:12:37 2023

@author: robbi
"""

import configparser
import requests
from icecream import ic
import ast

def get_key():
    config = configparser.ConfigParser()
    config.read('G:/My Drive/Data/data/open_ai/config.ini')  # Use the full path to the config file
    api_key = config['API']['key']
    return api_key


def get_response(question):
    openai_url = 'https://api.openai.com/v1/chat/completions'
    key = get_key()
    
    headers = {
        "Authorization": f"Bearer {key}",
        "Content-Type": "application/json"
    }

    data = {
        "model": "gpt-4",
        "messages": [{"role": "system", "content": "Your answeres are just lists containing the requested human uniprot proteins in the following format ['x','y']"},
                     {"role": "user", "content": question}],
        "max_tokens": 500  # Adjust as needed
    }
    
    try:
        response = requests.post(openai_url, json=data, headers=headers)
        if response.status_code == 200:
            response_data = response.json()
            if "choices" in response_data and len(response_data["choices"]) > 0:
                if "message" in response_data["choices"][0]:
                    response_str = response_data["choices"][0]["message"]["content"]
                    try:
                        
                        response_list = ast.literal_eval(response_str)
                        response = [s.upper() for s in response_list]
                        print(response)
                        return response
                    except ValueError:
                        print("Error in converting response to list.")
                else:
                    print("No 'message' key in the response.")
            else:
                print("No 'choices' in the response or 'choices' array is empty.")
        else:
            print("Error:", response.status_code, response.text)
    except Exception as e:
        print(f"An error occurred: {e}")
        
    # response = requests.post(openai_url, json=data, headers=headers)
    # if response.status_code == 200:
    #     response_data = response.json()
    #     if "choices" in response_data and len(response_data["choices"]) > 0:
    #         if "message" in response_data["choices"][0]:
    #             response_list = response_data["choices"][0]["message"]["content"]
    #             response_list = ast.literal_eval(response_list)
    #             return response_list
    #         else:
    #             print("No 'message' key in the response.")
    #     else:
    #         print("No 'choices' in the response or 'choices' array is empty.")
    # else:
    #     print("Error:", response.status_code, response.text)

def get_list(question):
    proteins = get_response(question)
    return proteins
        
question = "Give me a list of FBXO Family (F-box only proteins)"

proteins = get_response(question)
