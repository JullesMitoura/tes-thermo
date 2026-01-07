"""
ThermoAgent using OpenAI SDK directly (without LangChain).
Implements function calling natively with OpenAI API.
"""
from openai import AzureOpenAI
from typing import List, Dict, Any
from tes_thermo.thermo_agent.ming_function import min_g_function, get_min_g_function_schema
from tes_thermo.utils.prompts import Prompts


class Agent:
    """
    Agent using OpenAI SDK directly for function calling.
    Replaces LangChain AgentExecutor with native OpenAI implementation.
    """
    
    def __init__(self, 
                 openai_client,
                 llm_model_name: str,
                 embedding_model_name: str = None,
                 vsearch=None):
        """
        Initialize Agent.
        
        Args:
            openai_client: OpenAI client instance (OpenAI or AzureOpenAI)
            llm_model_name: Name of the LLM model to use
            embedding_model_name: Name of the embedding model (optional, for RAG)
            vsearch: VectorSearch instance (optional, for RAG)
        """
        self.llm = openai_client
        self.llm_model_name = llm_model_name
        self.embedding_model_name = embedding_model_name
        self.vsearch = vsearch
        self.system_prompt = Prompts.thermo_agent()
        
        # Define available functions
        self.functions = {
            "min_g_calc": min_g_function
        }
        
        # Function schemas for OpenAI
        self.function_schemas = [
            get_min_g_function_schema()
        ]
    
    def _perform_rag_search(self, query: str) -> str:
        """Perform RAG search and return context as string."""
        if self.vsearch is None or self.embedding_model_name is None:
            return ""
        try:
            docs = self.vsearch.search(query=query, k=10)
            if docs:
                context = " ".join(doc.get("text", "") for doc in docs)
                return f"\n\n[Context from documents:]\n{context}\n"
            return ""
        except Exception as e:
            print(f"Error during RAG search: {e}")
            return ""
    
    def _call_function(self, function_name: str, arguments: Dict[str, Any]) -> str:
        """
        Call a function by name with arguments.
        
        Args:
            function_name: Name of the function to call
            arguments: Arguments to pass to the function
        
        Returns:
            Function result as string
        """
        if function_name in self.functions:
            try:
                result = self.functions[function_name](**arguments)
                return str(result)
            except Exception as e:
                return f"Error calling function {function_name}: {str(e)}"
        else:
            return f"Function {function_name} not found"
    
    def run(self, messages: List[Dict[str, str]]) -> Dict[str, Any]:
        """
        Run the agent with a conversation.
        
        Args:
            messages: List of message dictionaries with 'role' and 'content' keys
        
        Returns:
            Dictionary with 'output' key containing the final response
        """
        # Add RAG context to the last user message if available
        if self.vsearch and messages:
            last_message = messages[-1]
            if last_message.get("role") == "user":
                rag_context = self._perform_rag_search(last_message.get("content", ""))
                if rag_context:
                    last_message["content"] = rag_context + "\n\n" + last_message["content"]
        
        # Prepare messages with system prompt
        conversation_messages = [
            {"role": "system", "content": self.system_prompt}
        ] + messages
        
        max_iterations = 10
        iteration = 0
        
        while iteration < max_iterations:
            iteration += 1
            
            try:
                # Call OpenAI API
                # Extract function schemas (OpenAI expects "type": "function" format)
                tools = [
                    {
                        "type": "function",
                        "function": schema["function"]
                    }
                    for schema in self.function_schemas
                ]
                
                response = self.llm.chat.completions.create(
                    model=self.llm_model_name,
                    messages=conversation_messages,
                    tools=tools if tools else None,
                    tool_choice="auto" if tools else None
                )
                
                message = response.choices[0].message
                
                # Check if function calling is required
                tool_calls = message.tool_calls
                
                if not tool_calls:
                    # No function calls, return the response
                    return {
                        "output": message.content or ""
                    }
                
                # Add assistant message with tool_calls to conversation
                assistant_message = {
                    "role": "assistant",
                    "content": message.content or None
                }
                
                # Include tool_calls in assistant message
                if tool_calls:
                    assistant_message["tool_calls"] = [
                        {
                            "id": tool_call.id,
                            "type": "function",
                            "function": {
                                "name": tool_call.function.name,
                                "arguments": tool_call.function.arguments
                            }
                        }
                        for tool_call in tool_calls
                    ]
                
                conversation_messages.append(assistant_message)
                
                # Process tool calls
                for tool_call in tool_calls:
                    function_name = tool_call.function.name
                    try:
                        import json
                        arguments = json.loads(tool_call.function.arguments)
                    except json.JSONDecodeError:
                        arguments = {}
                    
                    # Call the function
                    function_result = self._call_function(function_name, arguments)
                    
                    # Add function result to conversation (OpenAI expects "tool" role with tool_call_id)
                    conversation_messages.append({
                        "role": "tool",
                        "tool_call_id": tool_call.id,
                        "content": function_result
                    })
                
                # Continue loop to get final response
                
            except Exception as e:
                return {
                    "output": f"Error during agent execution: {str(e)}"
                }
        
        # If we've reached max iterations, return the last message
        if conversation_messages:
            last_msg = conversation_messages[-1]
            return {
                "output": last_msg.get("content", "Maximum iterations reached")
            }
        
        return {
            "output": "No response generated"
        }
