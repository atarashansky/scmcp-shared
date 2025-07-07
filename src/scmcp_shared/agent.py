import instructor
from openai import OpenAI
from scmcp_shared.schema.tool import ToolList
import os


from agno.agent import Agent
from agno.models.openai import OpenAILike
from scmcp_shared.kb import load_kb

model = OpenAILike(
    id=os.getenv("MODEL"),
    base_url=os.getenv("BASE_URL"),
    api_key=os.getenv("API_KEY"),
)


def rag_agent(task, software=None):
    knowledge_base = load_kb(software=software)
    agent = Agent(
        model=model,
        knowledge=knowledge_base,
        show_tool_calls=True,
        search_knowledge=True,
    )
    query = f"""
    <task>
    {task}
    </task>
    查询知识库，给出一个用于解决任务的代码示例。返回结果格式为：
    <code_example>
        [code_example]
    </code_example>
    """
    rep = agent.run(query)
    return rep.content


def select_tool(query):
    API_KEY = os.environ.get("API_KEY", None)
    BASE_URL = os.environ.get("BASE_URL", None)
    MODEL = os.environ.get("MODEL", None)

    client = OpenAI(api_key=API_KEY, base_url=BASE_URL)
    client = instructor.from_openai(client)

    response = client.chat.completions.create(
        model=MODEL,
        messages=[
            {
                "role": "system",
                "content": "you are a bioinformatician, you are given a task and a list of tools, you need to select the most directly relevant tools to use to solve the task",
            },
            {"role": "user", "content": query},
        ],
        response_model=ToolList,
    )
    return response.tools
