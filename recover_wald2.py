import json

with open('extract.jsonl', 'r') as f:
    for line in f:
        data = json.loads(line)
        content = str(data)
        if 'struct split_lognormal_shape' in content:
            print("Found in step type:", data.get('type'), "source:", data.get('source'))
            # We can extract the content string if it's in a specific field
            if 'tool_responses' in data:
                for resp in data['tool_responses']:
                    if 'struct split_lognormal_shape' in str(resp):
                        with open('recovered.txt', 'w') as out:
                            out.write(str(resp))
            elif 'content' in data:
                with open('recovered.txt', 'w') as out:
                    out.write(data['content'])
