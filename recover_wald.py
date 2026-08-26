import json

with open('extract.jsonl', 'r') as f:
    lines = f.readlines()

for line in lines:
    try:
        data = json.loads(line)
        # Look for a tool response for view_file
        if 'tool_calls' in data:
            # wait, tool outputs are in another step type or in content?
            pass
        if data.get('type') == 'SYSTEM_RESPONSE' or data.get('type') == 'PLANNER_RESPONSE' or 'content' in data:
            content = data.get('content', '')
            if 'split_lognormal_shape' in content and 'wald_functions.h' in content:
                # Let's just find the file content.
                # Usually view_file returns standard output format.
                pass
    except Exception as e:
        pass

# Let's just parse the full transcript to find the exact view_file output.
with open('/home/ubuntu/.gemini/antigravity-cli/brain/aac58f35-665d-4482-af2c-ba17375b921a/.system_generated/logs/transcript_full.jsonl', 'r') as f:
    for line in f:
        data = json.loads(line)
        if data.get('type') == 'ACTION_RESULT':
            content = data.get('content', '')
            if 'split_lognormal_shape_params' in content and 'struct split_lognormal_shape' in content:
                with open('recovered_wald_functions.h.txt', 'w') as out:
                    out.write(content)
                print("Recovered to recovered_wald_functions.h.txt")
                break
