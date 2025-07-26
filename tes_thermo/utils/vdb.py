import faiss
import fitz
from langchain.text_splitter import MarkdownTextSplitter
from langchain_community.docstore.in_memory import InMemoryDocstore
from langchain_community.vectorstores import FAISS
from langchain_core.documents import Document
from langchain_openai import AzureOpenAIEmbeddings
from typing import Dict, List, Optional
import os


class DocumentProcessor:
    """
    Helper class responsible for extracting and processing text from documents.
    Provides text filtering, PDF extraction, and chunking for embedding preparation.
    """
    
    def _filter_extracted_text(self, text: str) -> str:
        """Filter out empty lines or lines with only one word."""
        if not text:
            return ""
        
        lines = text.splitlines()
        filtered_lines = [
            line for line in lines
            if line.strip() and len(line.strip()) > 1 and len(line.strip().split()) > 1
        ]
        return "\n".join(filtered_lines)
    
    def extract_pdf(self, file_bytes: bytes) -> str:
        """Extract text from PDF bytes and filter the output."""
        try:
            doc = fitz.open(stream=file_bytes, filetype="pdf")
            full_text = "".join(page.get_text("text") for page in doc)
            doc.close()
            return self._filter_extracted_text(full_text)
        except Exception as e:
            print(f"Error while processing PDF from bytes: {e}")
            return ""
    
    def create_chunks(self, text: str, chunk_size: int = 1400, chunk_overlap: int = 200) -> List[str]:
        """Split text into chunks using a MarkdownTextSplitter."""
        if not text:
            return []
        
        markdown_splitter = MarkdownTextSplitter(
            chunk_size=chunk_size, 
            chunk_overlap=chunk_overlap
        )
        return markdown_splitter.split_text(text)


class VectorSearch:
    """
    Wrapper class for FAISS vector search integration with Azure OpenAI embeddings.
    Handles the creation of the FAISS index and performs similarity search.
    """
    
    def __init__(self, vector_store: FAISS):
        self.vector_store = vector_store
    
    @classmethod
    def from_documents(cls,
                      document_paths: List[str],
                      embedding: AzureOpenAIEmbeddings,
                      dimension: int = 1536) -> "VectorSearch":
        """
        Build a FAISS vector store from a list of document file paths.
        
        Args:
            document_paths (List[str]): List of file paths to PDF documents.
            embedding (AzureOpenAIEmbeddings): Embedding model instance.
            dimension (int): Dimension of the embeddings.
        
        Returns:
            VectorSearch: Instance with the created FAISS store.
        """
        processor = DocumentProcessor()
        all_chunks_for_db = {}
        
        # Process each document path
        for path in document_paths:
            try:
                # Use filename as description
                description = os.path.basename(path)
                
                # Read file content as bytes
                with open(path, "rb") as f:
                    doc_bytes = f.read()
                    
            except FileNotFoundError:
                print(f"File not found: '{path}'. Skipping.")
                continue
            except Exception as e:
                print(f"Error reading file '{path}': {e}")
                continue
            
            # Extract text from PDF
            text = processor.extract_pdf(doc_bytes)
            if not text:
                print(f"⚠️ No text extracted from '{path}'. Skipping.")
                continue
            
            # Create chunks from extracted text
            chunks = processor.create_chunks(text)
            
            # Generate source ID from filename
            source_id = description.replace(' ', '_').lower().removesuffix('.pdf')
            
            # Add chunks to database dictionary
            for i, chunk in enumerate(chunks):
                chunk_id = f"{source_id}_{i}"
                all_chunks_for_db[chunk_id] = {
                    "text": chunk, 
                    "source": description
                }
        
        # Create FAISS index
        index = faiss.IndexFlatL2(dimension)
        vector_store = FAISS(
            embedding_function=embedding,
            index=index,
            docstore=InMemoryDocstore(),
            index_to_docstore_id={},
        )
        
        # Prepare documents for vector store
        documents = []
        ids = []
        
        for chunk_id, chunk_data in all_chunks_for_db.items():
            doc = Document(
                page_content=chunk_data["text"],
                metadata={"source": chunk_data.get("source", "unknown")}
            )
            documents.append(doc)
            ids.append(chunk_id)
        
        # Add documents to vector store
        if documents:
            vector_store.add_documents(documents=documents, ids=ids)
            print(f"Index created successfully with {len(documents)} chunks.")
        else:
            print("No documents were processed. Empty index created.")
        
        return cls(vector_store)
    
    def search(self,
               query: str,
               k: int = 10,
               filter: Optional[str] = None) -> List[Document]:
        """
        Perform a similarity search in the FAISS index.
        
        Args:
            query (str): Search query text.
            k (int): Number of results to return.
            filter (Optional[str]): Optional filter by document source.
        
        Returns:
            List[Document]: List of retrieved documents.
        """
        if not self.vector_store:
            return []
        
        search_kwargs = {'k': k}
        if filter:
            search_kwargs['filter'] = {"source": filter}
        
        try:
            return self.vector_store.similarity_search(query, **search_kwargs)
        except Exception as e:
            print(f"Error during search: {e}")
            return []