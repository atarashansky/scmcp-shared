import json
from fastmcp import FastMCP
from fastmcp.server.dependencies import get_context
from ..agent import select_tool
from pydantic import Field
from typing import Union

auto_mcp = FastMCP("SmartMCP-select-Server")


@auto_mcp.tool(tags={"auto"})
def search_tool(
    task: str = Field(
        description="The tasks or questions that needs to be solved using available tools"
    ),
):
    """search the tools and get tool parameters that can be used to solve the  user's tasks or questions"""
    ctx = get_context()
    ads = ctx.request_context.lifespan_context
    adata = ads.get_adata()
    fastmcp = ctx.fastmcp
    if hasattr(fastmcp._tool_manager, "_all_tools"):
        all_tools = fastmcp._tool_manager._all_tools
    else:
        all_tools = fastmcp._tool_manager._tools
    auto_tools = fastmcp._tool_manager._tools
    fastmcp._tool_manager._tools = all_tools
    query = f"""
    <adata>
        {str(adata)}
    </adata>
    <task>{task}</task>\n
    """
    for name in all_tools:
        tool = all_tools[name]
        query += f"<Tool>\n<name>{name}</name>\n<description>{tool.description}</description>\n</Tool>\n"
    results = select_tool(query)
    tool_list = []
    for tool in results.tools:
        tool = tool.model_dump()
        tool["parameters"] = all_tools[tool["name"]].parameters
        tool_list.append(tool)
    fastmcp._tool_manager._tools = auto_tools
    return tool_list


@auto_mcp.tool(tags={"auto"})
async def run_tool(
    name: str = Field(description="The name of the tool to run"),
    parameter: Union[dict, str] = Field(description="The parameters to pass to the tool"),
):
    """run the tool with the given name and parameters. Only start call the tool when last tool is finished."""
    ctx = get_context()
    fastmcp = ctx.fastmcp
    all_tools = fastmcp._tool_manager._all_tools
    auto_tools = fastmcp._tool_manager._tools
    fastmcp._tool_manager._tools = all_tools

    # Deserialize JSON string if needed
    if isinstance(parameter, str):
        try:
            parameter = json.loads(parameter)
        except json.JSONDecodeError:
            fastmcp._tool_manager._tools = auto_tools
            return {"error": f"Invalid JSON in parameter: {parameter}"}

    try:
        result = await fastmcp._tool_manager.call_tool(name, parameter)
    except Exception as e:
        fastmcp._tool_manager._tools = auto_tools
        result = {"error": str(e)}
    fastmcp._tool_manager._tools = auto_tools
    return result
