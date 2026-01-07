"""
ThermoAgent main interface.
Uses OpenAI SDK directly (no LangChain).
"""
from tes_thermo.thermo_agent.agent import Agent
from typing import Dict, Any


class ThermoAgent:
    """
    Main interface for ThermoAgent.
    Simplified to work with OpenAI SDK directly.
    """
    
    def __init__(self, 
                 openai_client,
                 llm_model_name: str,
                 embedding_model_name: str = None,
                 vsearch=None):
        """
        Initialize ThermoAgent.
        
        Args:
            openai_client: OpenAI client instance (OpenAI or AzureOpenAI)
            llm_model_name: Name of the LLM model to use
            embedding_model_name: Name of the embedding model (optional, for RAG)
            vsearch: VectorSearch instance (optional, for RAG)
        """
        self.vsearch = vsearch
        self.agent = Agent(
            openai_client=openai_client,
            llm_model_name=llm_model_name,
            embedding_model_name=embedding_model_name,
            vsearch=self.vsearch
        )
        self.chat_history = []
    
    def chat(self, prompt: str) -> Dict[str, Any]:
        """
        Chat with the agent.
        
        Args:
            prompt: User prompt string
        
        Returns:
            Dictionary with 'output' key containing the response
        """
        # Add user message to history
        self.chat_history.append({
            "role": "user",
            "content": prompt
        })
        
        # Run agent
        result = self.agent.run(messages=self.chat_history)
        
        # Add assistant response to history
        if "output" in result:
            self.chat_history.append({
                "role": "assistant",
                "content": result["output"]
            })
        
        return result
