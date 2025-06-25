import os
import time

import openai
from dotenv import load_dotenv

load_dotenv()

OPENAI_API_KEY = os.environ.get("OPENAI_API_KEY")
ASSISTANT_ID = os.environ.get("OPENAI_ASSISTANT_ID")

client = openai.OpenAI(api_key=OPENAI_API_KEY)
thread = client.beta.threads.create()

message = client.beta.threads.messages.create(
    thread_id=thread.id,
    role="user",
    content="53 years old, 90 kg, allergic to penicillin, undergoing implant surgery.",
)


run = client.beta.threads.runs.create(thread_id=thread.id, assistant_id=ASSISTANT_ID)

while True:
    run_status = client.beta.threads.runs.retrieve(thread_id=thread.id, run_id=run.id)
    if run_status.status == "completed":
        break
    time.sleep(1)

messages = client.beta.threads.messages.list(thread_id=thread.id)
print(messages.data[0].content[0].text.value)
