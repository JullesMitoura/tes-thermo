from langchain.tools import BaseTool
from tes_thermo.utils.prompts import Prompts
from tes_thermo.utils import VectorSearch
from langchain_openai import AzureOpenAIEmbeddings
from typing import Optional
from pydantic import BaseModel, Field

class RaGGInputs(BaseModel):
    query: Optional[str] = Field(description="Based on the user's question, define a query to be made in the vector database.")

class RaG(BaseTool):
    name: str = "rag"
    description: str = Prompts.rag()
    args_schema = RaGGInputs
    embedding: AzureOpenAIEmbeddings
    vsearch: VectorSearch

    def _run(self, 
             query: str) -> str:
        docs = self.vsearch.search(query = query)
        text = " ".join(doc.page_content for doc in docs)
        return text