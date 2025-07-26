from langchain_core.prompts import ChatPromptTemplate, MessagesPlaceholder
from langchain.agents import AgentExecutor, create_openai_tools_agent

class CreateAgent():
    def __init__(self,
                 llm,
                 tools,
                 system_prompt):
        
        self.llm = llm
        self.tools = tools
        self.system_prompt = system_prompt

    def create_node(self):
        prompt = ChatPromptTemplate.from_messages(
            [
                (
                    "system",
                    self.system_prompt,
                ),
                MessagesPlaceholder(variable_name="messages"),
                MessagesPlaceholder(variable_name="agent_scratchpad"),
            ]
        )
        agent = create_openai_tools_agent(self.llm, self.tools, prompt)
        executor = AgentExecutor(agent=agent, 
                                 tools=self.tools)
        return executor