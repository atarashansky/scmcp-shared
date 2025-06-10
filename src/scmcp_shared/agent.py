import instructor
from openai import OpenAI
from scmcp_shared.schema.tool import ToolList
import os


client = OpenAI(api_key=os.environ.get("API_KEY"), base_url=os.environ.get("BASE_URL"))
client = instructor.from_openai(client)


def select_tool(query):

    response = client.chat.completions.create(
        model=os.environ.get("MODEL"),
        messages=[
            {
                "role": "system", 
                "content": f"you are a bioinformatician, you are given a task and a list of tools, you need to select the most directly relevant tools to use to solve the task"},
            {
                "role": "user",
                "content": query
            },  
        ],
        response_model=ToolList,
    )
    return response.tools
