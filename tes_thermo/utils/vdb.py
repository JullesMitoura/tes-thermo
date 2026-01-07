"""
Vector database using FAISS nativo (sem LangChain).
"""
import faiss
import fitz
import numpy as np
from typing import List, Optional, Dict, Any
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
        """
        Split text into chunks using simple text splitting.
        Replaces LangChain MarkdownTextSplitter with native implementation.
        """
        if not text:
            return []
        
        # Simple chunking implementation
        words = text.split()
        chunks = []
        
        if len(words) <= chunk_size:
            return [text]
        
        # Split into chunks with overlap
        start = 0
        while start < len(words):
            end = min(start + chunk_size, len(words))
            chunk = " ".join(words[start:end])
            chunks.append(chunk)
            
            if end >= len(words):
                break
                
            # Move start by chunk_size - overlap
            start += chunk_size - chunk_overlap
        
        return chunks


class VectorSearch:
    """
    Vector search using FAISS nativo without LangChain.
    Uses OpenAI embeddings directly.
    """
    
    def __init__(self, 
                 index: faiss.Index,
                 documents: List[Dict[str, Any]],
                 embeddings: np.ndarray,
                 embedding_function):
        """
        Initialize VectorSearch with FAISS index and documents.
        
        Args:
            index: FAISS index
            documents: List of document dictionaries with 'text' and 'metadata' keys
            embeddings: Numpy array of embeddings
            embedding_function: Function to generate embeddings
        """
        self.index = index
        self.documents = documents
        self.embeddings = embeddings
        self.embedding_function = embedding_function
    
    @classmethod
    def from_documents(cls,
                      document_paths: List[str],
                      openai_client,
                      model_name: str,
                      dimension: int = 1536) -> "VectorSearch":
        """
        Build a FAISS vector store from a list of document file paths.
        
        Args:
            document_paths: List of file paths to PDF documents
            openai_client: OpenAI client instance (OpenAI or AzureOpenAI)
            model_name: Name of the embedding model to use
            dimension: Dimension of the embeddings
        
        Returns:
            VectorSearch: Instance with the created FAISS store
        """
        processor = DocumentProcessor()
        all_chunks = []
        
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
            
            # Add chunks to list
            for chunk in chunks:
                all_chunks.append({
                    "text": chunk,
                    "metadata": {"source": description}
                })
        
        if not all_chunks:
            print("No documents were processed. Creating empty index.")
            # Create empty index
            index = faiss.IndexFlatL2(dimension)
            
            def embedding_function(texts):
                if isinstance(texts, str):
                    texts = [texts]
                response = openai_client.embeddings.create(
                    model=model_name,
                    input=texts
                )
                return [item.embedding for item in response.data]
            
            return cls(index, [], np.array([]), embedding_function)
        
        # Generate embeddings for all chunks
        print(f"Generating embeddings for {len(all_chunks)} chunks...")
        texts = [chunk["text"] for chunk in all_chunks]
        
        # Create embedding function
        def embedding_function(texts):
            if isinstance(texts, str):
                texts = [texts]
            response = openai_client.embeddings.create(
                model=model_name,
                input=texts
            )
            return [item.embedding for item in response.data]
        
        embeddings = embedding_function(texts)
        
        # Convert to numpy array
        embeddings_array = np.array(embeddings, dtype=np.float32)
        
        # Create FAISS index
        index = faiss.IndexFlatL2(dimension)
        index.add(embeddings_array)
        
        print(f"Index created successfully with {len(all_chunks)} chunks.")
        
        return cls(index, all_chunks, embeddings_array, embedding_function)
    
    def search(self,
               query: str,
               k: int = 10,
               filter: Optional[str] = None) -> List[Dict[str, Any]]:
        """
        Perform a similarity search in the FAISS index.
        
        Args:
            query: Search query text
            k: Number of results to return
            filter: Optional filter by document source
        
        Returns:
            List of document dictionaries with 'text' and 'metadata' keys
        """
        if not self.index or self.index.ntotal == 0:
            return []
        
        try:
            # Generate embedding for query
            query_embedding = self.embedding_function([query])
            query_vector = np.array(query_embedding, dtype=np.float32)
            
            # Search in FAISS index
            distances, indices = self.index.search(query_vector, k)
            
            # Get results
            results = []
            for i, idx in enumerate(indices[0]):
                if idx < len(self.documents):
                    doc = self.documents[idx]
                    # Apply filter if specified
                    if filter is None or doc["metadata"].get("source") == filter:
                        results.append({
                            "text": doc["text"],
                            "metadata": doc["metadata"],
                            "distance": float(distances[0][i])
                        })
            
            return results
            
        except Exception as e:
            print(f"Error during search: {e}")
            return []
