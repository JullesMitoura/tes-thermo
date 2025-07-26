from tes_thermo.thermo_agent.agent import Agent
from langchain_core.messages import HumanMessage, AIMessage

class ThermoAgent:
    def __init__(self, 
                 llm, 
                 embedding_model,
                 vsearch = None):
        self.vsearch = vsearch
        self.agent = Agent(llm=llm,
                           embedding_model = embedding_model,
                           vsearch = self.vsearch)
        self.chat_history = []

    def chat(self, prompt):
        self.chat_history.append(HumanMessage(content=prompt))
        result = self.agent.run(conversation_input = {"messages": self.chat_history})
        self.chat_history.append(AIMessage(content=result['output']))
        return result