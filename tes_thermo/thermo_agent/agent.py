from langchain_openai import AzureChatOpenAI
from langchain_core.messages import HumanMessage
from tes_thermo.thermo_agent.create_agents import CreateAgent
from tes_thermo.thermo_agent.ming_tool import MinG
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

        self.agent = CreateAgent(
            llm=self.llm,
            tools=tools,
            system_prompt=Prompts.thermo_agent()
        ).create_node()

    def _has_documents(self):
        """Check if vsearch has documents indexed."""
        if self.vsearch is None:
            return False
        try:
            # Try to get a sample search to see if there are documents
            # If the vector store is empty, it will return empty list
            results = self.vsearch.search("test", k=1)
            return len(results) > 0
        except:
            return False
    
    def _perform_rag_search(self, query: str) -> str:
        """Perform RAG search and return context as string."""
        if self.vsearch is None:
            return ""
        try:
            docs = self.vsearch.search(query=query, k=10)
            if docs:
                context = " ".join(doc.page_content for doc in docs)
                return f"\n\n[Context from documents:]\n{context}\n"
            return ""
        except Exception as e:
            print(f"Error during RAG search: {e}")
            return ""

    def run(self, conversation_input: dict):
        # Perform automatic RAG if documents are available
        if self._has_documents():
            # Get the last user message (HumanMessage)
            messages = conversation_input.get("messages", [])
            user_query = None
            last_user_msg_idx = None
            
            # Find the last HumanMessage
            for idx, msg in enumerate(messages):
                if isinstance(msg, HumanMessage) and hasattr(msg, 'content') and msg.content:
                    user_query = msg.content
                    last_user_msg_idx = idx
            
            if user_query and last_user_msg_idx is not None:
                rag_context = self._perform_rag_search(user_query)
                if rag_context:
                    # Prepend RAG context to the last user message
                    last_msg = messages[last_user_msg_idx]
                    last_msg.content = rag_context + "\n\n" + last_msg.content
        
        return self.agent.invoke(conversation_input)