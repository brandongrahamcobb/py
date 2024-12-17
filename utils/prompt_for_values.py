def prompt_for_values(prompt: str, default_value: str) -> str:
    value = input(f'{prompt} [{default_value}]: ')
    return value if value else default_value
