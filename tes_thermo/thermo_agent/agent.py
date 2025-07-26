from langchain_openai import AzureChatOpenAI
from tes_thermo.thermo_agent.create_agents import CreateAgent
from tes_thermo.thermo_agent.tools.ming_tool import MinG
from tes_thermo.thermo_agent.tools.search_tool import RaG
from tes_thermo.utils.prompts import Prompts

class Agent:
    def __init__(self, 
                 llm: AzureChatOpenAI,
                 embedding_model,
                 vsearch=None):

        self.llm = llm
        self.embedding_model = embedding_model
        self.vsearch = vsearch

        tools = [MinG()]
        if self.vsearch is not None:
            tools.append(RaG(vsearch=self.vsearch, embedding=self.embedding_model))

        self.agent = CreateAgent(
            llm=self.llm,
            tools=tools,
            system_prompt=Prompts.thermo_agent()
        ).create_node()

    def run(self, conversation_input: dict):
        return self.agent.invoke(conversation_input)