# Локальный офлайн-агент на free-моделях: пошаговая инструкция

Ниже — рабочая архитектура и код для агента, который:

1. Работает полностью офлайн (модель, память, скиллы, RAG — всё на диске);
2. Вызывается из консоли (`agent ...`);
3. Умеет **создавать** новые скиллы (пишет и регистрирует Python-модули);
4. Умеет **пользоваться** скиллами (вызывает их как tool calls);
5. При явном указании документа ищет **только в нём**, а не по всей базе;
6. Выходит в интернет только через браузер, поднятый **в отдельном Docker-контейнере**;
7. Пишет и исполняет код (Python, C++) в песочнице, с опциональной работой из VSCode.

## Архитектура

```
┌───────────────────────────────────────────────────────────────────────────┐
│                             agent (CLI)                                   │
│                Python-процесс: REPL / однократный вызов                   │
├──────────────┬──────────────┬──────────────┬──────────────┬───────────────┤
│  Ollama      │  Skills      │  RAG         │  Sandbox     │  Memory       │
│  (LLM        │  registry    │  (Chroma +   │  (Docker:    │  (sessions/,  │
│  runtime,    │  + create_   │  приоритет   │  python,     │  суммаризация │
│  localhost:  │  skill       │  документа)  │  gcc,        │  каждые 5 пар │
│  11434)      │              │              │  browser)    │  реплик)      │
└──────────────┴──────────────┴──────────────┴──────────────┴───────────────┘
```

Все компоненты, кроме контейнера-браузера, никогда не открывают сетевые соединения наружу — это и даёт офлайн-режим "по умолчанию". Блок Memory физически — папка `sessions/` на диске плюс код в `agent/memory.py` (Шаг 8); он появился позже остальных при доработке инструкции, поэтому упоминается отдельно, а не встроен в исходную схему из четырёх блоков.

## 0. Требования к железу

- 16 ГБ RAM минимум (32 ГБ комфортно), GPU с 8–12 ГБ VRAM ускоряет генерацию, но не обязателен — на CPU модели 7B–14B в квантовании Q4 работают, просто медленнее.
- Docker Desktop (или Docker Engine на Linux) — для песочницы кода и браузера.
- Python 3.11+.
- ~15–40 ГБ диска под модели.

## Шаг 1. Установить Ollama и модель(и)

Ollama — локальный runtime для open-моделей, отдаёт OpenAI-совместимый API на `http://localhost:11434`, после скачивания модели не требует интернета.

```bash
# macOS / Linux
curl -fsSL https://ollama.com/install.sh | sh
# Windows: скачать установщик с ollama.com

ollama serve            # поднять сервер (на Linux обычно уже сервис)
```

Скачайте модели заранее (это единственный момент, когда нужен интернет):

```bash
# основная модель для рассуждений и работы со скиллами (умеет tool calling)
ollama pull qwen2.5:14b-instruct-q4_K_M

# модель, заточенная под код (Python/C++)
ollama pull qwen2.5-coder:14b-instruct-q4_K_M

# embedding-модель для RAG
ollama pull nomic-embed-text
```

Если видеокарты нет или слабая — берите `:7b` вместо `:14b`. Проверить, что всё работает офлайн: выключите Wi‑Fi/Ethernet и выполните `ollama run qwen2.5:14b-instruct-q4_K_M "привет"` — ответ должен прийти без интернета.

## Шаг 2. Каркас проекта

```bash
mkdir -p ~/local-agent/{skills,docs,sessions,tmp,sandbox_images}
cd ~/local-agent
python3 -m venv .venv && source .venv/bin/activate
pip install ollama chromadb pyyaml docker typer rich
```

Структура:

```
local-agent/
├── agent/
│   ├── main.py          # REPL / CLI
│   ├── core.py          # цикл рассуждений + вызов моделей
│   ├── skills.py        # реестр скиллов + create_skill
│   ├── skill_guard.py   # проверка кода скиллов перед записью
│   ├── memory.py        # сохранение и суммаризация сессий
│   ├── rag.py           # индексация и поиск по документам
│   └── sandbox.py       # запуск кода и браузера в Docker
├── skills/              # скиллы (создаются агентом или вручную)
├── docs/                # документы для RAG — только источники для индексации
├── sessions/            # история диалогов (см. Шаг 8) — намеренно отдельно от docs/
└── tmp/                 # рабочие файлы песочницы
```

`docs/` и `sessions/` разведены умышленно: `docs/` — это источники, которые `ingest()` разбивает на чанки и кладёт в Chroma для RAG, а `sessions/` — это лог диалогов. Если сложить их в одну папку, есть риск, что случайный `ingest()` по всей `docs/` подхватит и файлы сессий, и тогда агент начнёт "вспоминать" куски старых разговоров как будто это содержимое учебника — то есть напрямую сломает требование про приоритет документа из Шага 4.

## Шаг 3. Реестр скиллов и авто-создание скиллов

Каждый скилл — папка с `manifest.yaml` (имя, описание, JSON-схема параметров) и `skill.py` с функцией `run(**kwargs) -> str`. Это тот же принцип, что у "skills" в современных агентных системах: манифест описывает, когда скилл нужен, а код — что он делает.

`agent/skills.py`:

```python
import importlib.util, yaml, pathlib, json

SKILLS_DIR = pathlib.Path(__file__).parent.parent / "skills"

def load_skills():
    skills = {}
    for folder in SKILLS_DIR.iterdir():
        manifest_path = folder / "manifest.yaml"
        if not manifest_path.exists():
            continue
        manifest = yaml.safe_load(manifest_path.read_text())
        spec = importlib.util.spec_from_file_location(
            manifest["name"], folder / "skill.py")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        skills[manifest["name"]] = {"manifest": manifest, "module": module}
    return skills

def to_ollama_tools(skills: dict) -> list:
    """Формат tool-описаний для передачи модели."""
    tools = []
    for name, s in skills.items():
        m = s["manifest"]
        tools.append({
            "type": "function",
            "function": {
                "name": name,
                "description": m["description"],
                "parameters": m.get("parameters", {"type": "object", "properties": {}}),
            },
        })
    return tools

def run_skill(skills: dict, name: str, args: dict) -> str:
    return skills[name]["module"].run(**args)


CREATE_SKILL_MANIFEST = {
    "name": "create_skill",
    "description": (
        "Создать новый скилл: сохраняет manifest.yaml и skill.py в папку skills/, "
        "делая новую возможность доступной для последующих вызовов."
    ),
    "parameters": {
        "type": "object",
        "properties": {
            "skill_name": {"type": "string"},
            "description": {"type": "string"},
            "parameters_schema": {"type": "object"},
            "python_code": {"type": "string", "description": "Тело файла skill.py с функцией run(**kwargs)"},
        },
        "required": ["skill_name", "description", "python_code"],
    },
}

def create_skill(skill_name, description, python_code, parameters_schema=None, force=False):
    from .skill_guard import validate_skill_code
    if not force:
        validate_skill_code(python_code)  # бросит SkillCodeError, если что-то запрещено

    folder = SKILLS_DIR / skill_name
    folder.mkdir(parents=True, exist_ok=True)
    manifest = {
        "name": skill_name,
        "description": description,
        "parameters": parameters_schema or {"type": "object", "properties": {}},
    }
    (folder / "manifest.yaml").write_text(yaml.dump(manifest, allow_unicode=True))
    (folder / "skill.py").write_text(python_code)
    return f"Скилл '{skill_name}' создан."
```

Проверка кода — `agent/skill_guard.py`. Важно не просто молча отклонять запрещённый код, а возвращать структурированную причину (какой именно импорт/вызов, в какой строке), чтобы модель могла объяснить это пользователю своими словами, а не сказать "не получилось":

```python
import ast

ALLOWED_IMPORTS = {
    "math", "datetime", "json", "re", "statistics", "textwrap",
    "collections", "itertools", "random", "string", "pathlib",
}
FORBIDDEN_CALLS = {"eval", "exec", "compile", "__import__"}
FORBIDDEN_ATTRS = {"system", "popen", "fork", "socket", "connect"}

class SkillCodeError(Exception):
    def __init__(self, restriction: str, detail: str):
        self.restriction = restriction  # машиночитаемый код ограничения, напр. "import:requests"
        self.detail = detail            # человекочитаемое объяснение для пользователя
        super().__init__(detail)

def validate_skill_code(code: str):
    tree = ast.parse(code)
    for node in ast.walk(tree):
        if isinstance(node, (ast.Import, ast.ImportFrom)):
            names = [n.name.split(".")[0] for n in node.names] if isinstance(node, ast.Import) \
                    else [node.module.split(".")[0]] if node.module else []
            for name in names:
                if name not in ALLOWED_IMPORTS:
                    raise SkillCodeError(
                        restriction=f"import:{name}",
                        detail=f"импорт модуля '{name}' запрещён для скиллов (строка {node.lineno})",
                    )
        if isinstance(node, ast.Call):
            fname = node.func.id if isinstance(node.func, ast.Name) else getattr(node.func, "attr", None)
            if fname in FORBIDDEN_CALLS or fname in FORBIDDEN_ATTRS:
                raise SkillCodeError(
                    restriction=f"call:{fname}",
                    detail=f"вызов '{fname}' запрещён для скиллов (строка {node.lineno})",
                )
    return True
```

`ALLOWED_IMPORTS` — это белый список, а не чёрный: любой модуль, которого там нет (в том числе `requests`, `os`, `subprocess`, `socket`), автоматически под запретом. Логика такая: если скиллу реально нужна сеть или процессы — для этого уже есть изолированные инструменты (`browse_web` в контейнере, `run_python`/`run_cpp` в песочнице), и создавать параллельный обходной путь через "обычный" скилл не стоит.

Дальше в Шаге 7 показано, как при срабатывании `SkillCodeError` агент не просто получает отказ, а формулирует пользователю, с каким конкретно ограничением он столкнулся, и ждёт решения — разрешить исключение или нет.

Пример готового скилла `skills/get_time/`:

`manifest.yaml`
```yaml
name: get_time
description: Возвращает текущие дату и время на машине пользователя.
parameters: {type: object, properties: {}}
```

`skill.py`
```python
from datetime import datetime

def run():
    return datetime.now().isoformat()
```

Агент вызывает `create_skill` как обычный tool, когда решает, что ему не хватает нужной способности — то есть "умеет составлять скиллы" реализуется как ещё один tool в его же цикле рассуждений, а не отдельная подсистема.

## Шаг 4. RAG с приоритетом явно указанного документа

`agent/rag.py`:

```python
import chromadb, ollama, pathlib

client = chromadb.PersistentClient(path="./chroma_db")
collection = client.get_or_create_collection("docs")

def embed(text: str):
    return ollama.embeddings(model="nomic-embed-text", prompt=text)["embedding"]

def ingest(path: str, chunk_size=800):
    text = pathlib.Path(path).read_text(errors="ignore")
    source = pathlib.Path(path).name
    chunks = [text[i:i+chunk_size] for i in range(0, len(text), chunk_size)]
    for idx, chunk in enumerate(chunks):
        collection.add(
            ids=[f"{source}-{idx}"],
            embeddings=[embed(chunk)],
            documents=[chunk],
            metadatas=[{"source": source}],
        )

def query(question: str, explicit_source: str | None = None, k=5):
    where = {"source": explicit_source} if explicit_source else None
    res = collection.query(query_embeddings=[embed(question)], n_results=k, where=where)
    return "\n---\n".join(res["documents"][0]) if res["documents"] else ""
```

Приоритет документа делаем простым и надёжным: не пытаемся угадывать источник по смыслу вопроса, а даём агенту явный параметр. В системном промпте перечисляем список загруженных источников (`collection.get()["metadatas"]`), и просим модель, если пользователь называет конкретный документ ("в учебнике Фихтенгольца", "по файлу lecture3.pdf"), — передавать `explicit_source` тем же именем файла в tool `search_docs`. Это снимает половину проблемы "модель домешивает данные из другого файла", а вторую половину закрывает `where`-фильтр в Chroma, который физически не даёт получить чанки из чужих документов.

`skills/search_docs/skill.py`:
```python
from agent.rag import query

def run(question: str, explicit_source: str = None):
    return query(question, explicit_source) or "Ничего не найдено."
```

## Шаг 5. Песочница для Python/C++

`agent/sandbox.py`:

```python
import subprocess, pathlib, textwrap, uuid

TMP = pathlib.Path("tmp")

def run_python(code: str) -> str:
    fid = uuid.uuid4().hex
    file = TMP / f"{fid}.py"
    file.write_text(code)
    result = subprocess.run([
        "docker", "run", "--rm", "--network", "none",
        "--memory", "512m", "--cpus", "1",
        "-v", f"{file.resolve()}:/sandbox/main.py:ro",
        "python:3.12-slim", "python", "/sandbox/main.py",
    ], capture_output=True, text=True, timeout=30)
    return result.stdout + result.stderr

def run_cpp(code: str) -> str:
    fid = uuid.uuid4().hex
    src = TMP / f"{fid}.cpp"
    src.write_text(code)
    cmd = f"g++ -O2 -o /sandbox/a.out /sandbox/{fid}.cpp && /sandbox/a.out"
    result = subprocess.run([
        "docker", "run", "--rm", "--network", "none",
        "--memory", "512m", "--cpus", "1",
        "-v", f"{TMP.resolve()}:/sandbox",
        "gcc:13", "bash", "-c", cmd,
    ], capture_output=True, text=True, timeout=30)
    return result.stdout + result.stderr
```

`--network none` физически отрезает контейнер от сети — код может ошибиться или быть недоверенным, но наружу он ничего не отправит. Оборачиваете это в скиллы `run_python` / `run_cpp` с соответствующими manifest.yaml, и агент получает возможность писать и тут же исполнять код.

## Шаг 6. Браузер в контейнере — единственный выход в интернет

Отдельный образ, у которого сеть есть, но он изолирован от остальной системы:

`sandbox_images/browser/Dockerfile`:
```dockerfile
FROM mcr.microsoft.com/playwright/python:v1.47.0
WORKDIR /app
COPY browse.py .
ENTRYPOINT ["python", "browse.py"]
```

`sandbox_images/browser/browse.py`:
```python
import sys
from playwright.sync_api import sync_playwright

url = sys.argv[1]
with sync_playwright() as p:
    browser = p.chromium.launch()
    page = browser.new_page()
    page.goto(url, timeout=20000)
    print(page.inner_text("body")[:5000])
    browser.close()
```

Сборка и вызов:
```bash
docker build -t agent-browser ./sandbox_images/browser
docker run --rm agent-browser "https://example.com"
```

Скилл `skills/browse_web/skill.py`:
```python
import subprocess

def run(url: str):
    r = subprocess.run(["docker", "run", "--rm", "agent-browser", url],
                        capture_output=True, text=True, timeout=40)
    return r.stdout or r.stderr
```

Важно: это единственный скилл, у контейнера которого не стоит `--network none`. Все остальные части агента (модель, RAG, скиллы, песочница кода) сеть не трогают — значит выключенный интернет не ломает ничего, кроме этого одного инструмента, и модель может сама сообщить пользователю, что сеть недоступна.

## Шаг 6.5. Вторая модель как скилл, а не отдельный собеседник

Ранее в инструкции скачивались две модели — `qwen2.5:14b-instruct` (общая, с tool calling) и `qwen2.5-coder:14b-instruct` (заточенная под код) — но в коде использовалась только первая. Это была недоработка: вторая модель нигде не вызывалась. Закрываю пробел через ещё один обычный скилл, `generate_code` — оркестратор (основная модель) решает вызвать его тем же способом, каким решает вызвать любой другой tool: по описанию в `manifest.yaml`.

`skills/generate_code/manifest.yaml`:
```yaml
name: generate_code
description: >
  Генерирует код на Python или C++ по техническому заданию с помощью модели,
  специализированной на коде. Обязательно используйте этот скилл для любой
  генерации кода длиннее нескольких строк, алгоритмов, структур данных или
  программ на C++ — не пишите код в ответе самостоятельно, а делегируйте сюда.
parameters:
  type: object
  properties:
    task: {type: string, description: "Техническое задание на естественном языке"}
    language: {type: string, enum: [python, cpp]}
  required: [task, language]
```

`skills/generate_code/skill.py`:
```python
import ollama

CODER_MODEL = "qwen2.5-coder:14b-instruct-q4_K_M"

def run(task: str, language: str = "python"):
    resp = ollama.chat(model=CODER_MODEL, messages=[
        {"role": "system", "content": (
            f"Ты — ассистент-программист. Пиши только код на {language}, "
            f"без пояснений вокруг и без markdown-разметки ```."
        )},
        {"role": "user", "content": task},
    ])
    return resp["message"]["content"]
```

Это делегирование, а не параллельный диалог: оркестратор вызывает `generate_code`, дожидается **одного** полного ответа coder-модели, кладёт его как обычный `tool`-результат обратно в свою историю сообщений — и продолжает цикл сам, при необходимости передавая полученный код дальше в `run_python`/`run_cpp` (Шаг 5). Coder-модель не видит остальной диалог с пользователем и не принимает решений — она чистая функция "задание → код", вызываемая по требованию.

Две модели по 14B, если обе окажутся загружены в память одновременно (Ollama выгружает неиспользуемые модели по таймауту, но при активном чередовании вызовов может держать обе), — это примерно вдвое больше RAM/VRAM, чем при одной. Если памяти не хватает, возьмите для coder-роли модель поменьше, например `qwen2.5-coder:7b-instruct-q4_K_M` — качество кода на уровне 7B-coder обычно выше, чем у 14B-generalist, так что даунгрейд именно этой модели относительно безопасен.

## Шаг 7. Основной цикл агента

Два важных момента, которых не было в первой версии: (1) если `create_skill` наткнулся на ограничение, решение "разрешить исключение" принимает **человек через терминал**, а не сама модель — иначе смысл whitelist теряется; (2) после каждого хода история проверяется на необходимость суммаризации (детали — в Шаге 8).

`agent/core.py`:

```python
import ollama, json
from . import skills as skills_mod
from .skill_guard import SkillCodeError
from .memory import maybe_summarize

MODEL = "qwen2.5:14b-instruct-q4_K_M"  # оркестратор; coder-модель вызывается отдельно внутри skills/generate_code

def confirm_restriction_with_user(error: SkillCodeError) -> bool:
    print(f"\n⚠️  Агент хочет создать скилл, но столкнулся с ограничением:")
    print(f"    {error.detail}")
    answer = input("    Разрешить это исключение для данного скилла? [y/N]: ").strip().lower()
    return answer == "y"

def run_agent(user_input: str, history: list):
    skills = skills_mod.load_skills()
    tools = skills_mod.to_ollama_tools(skills)
    tools.append({"type": "function", "function": skills_mod.CREATE_SKILL_MANIFEST})

    messages = history + [{"role": "user", "content": user_input}]
    while True:
        resp = ollama.chat(model=MODEL, messages=messages, tools=tools)
        msg = resp["message"]
        messages.append(msg)

        if not msg.get("tool_calls"):
            messages = maybe_summarize(messages, model=MODEL)
            return msg["content"], messages

        for call in msg["tool_calls"]:
            name = call["function"]["name"]
            args = call["function"]["arguments"]

            if name == "create_skill":
                try:
                    result = skills_mod.create_skill(**args)
                    skills = skills_mod.load_skills()
                except SkillCodeError as e:
                    if confirm_restriction_with_user(e):
                        result = skills_mod.create_skill(**args, force=True)
                        skills = skills_mod.load_skills()
                    else:
                        # Модель получает точную причину отказа и может
                        # объяснить её пользователю своими словами, а не
                        # просто сказать "не вышло".
                        result = (
                            f"Скилл не создан. Ограничение: {e.detail}. "
                            f"Пользователь не подтвердил исключение."
                        )
            else:
                result = skills_mod.run_skill(skills, name, args)

            messages.append({"role": "tool", "name": name, "content": str(result)})
```

`agent/main.py` (CLI) — с загрузкой и сохранением сессии (детали Шага 8):

```python
import typer, uuid, atexit
from .core import run_agent
from .memory import load_history, save_history

app = typer.Typer()
history = []
SESSION_ID = None

@app.command()
def chat(
    prompt: str = typer.Argument(None),
    session: str = typer.Option(None, "--session", help="Идентификатор сессии для сохранения/продолжения"),
    continue_last: bool = typer.Option(False, "--continue", help="Продолжить последнюю сессию"),
):
    global history, SESSION_ID
    SESSION_ID = session or ("last" if continue_last else uuid.uuid4().hex)
    history = load_history(SESSION_ID) if (session or continue_last) else []
    atexit.register(lambda: save_history(SESSION_ID, history))

    if prompt:
        answer, history = run_agent(prompt, history)
        print(answer)
    else:
        print(f"Локальный агент. Сессия: {SESSION_ID}. Ctrl+C для выхода.")
        while True:
            q = input("> ")
            answer, history = run_agent(q, history)
            print(answer)

if __name__ == "__main__":
    app()
```

`atexit.register(...)` гарантирует, что история будет сохранена на диск даже при закрытии консоли — не только по явному выходу из REPL, но и при Ctrl+C.

## Шаг 8. Память: суммаризация и хранение сессий

Первая версия хранила `history` только в оперативной памяти процесса, без ограничения размера, — при закрытии консоли всё терялось, а при долгом разговоре список рос бесконечно. Здесь это исправлено двумя механизмами.

**Где хранится сессия.** Логично положить её не в `docs/`, а в отдельную папку `sessions/` (см. обновлённую структуру в Шаге 2): `docs/` — это источники для RAG, которые `ingest()` без разбора режет на чанки и индексирует, а история диалога — совсем другая по смыслу и жизненному циклу сущность. Если держать их вместе, велик риск, что лог сессии случайно попадёт в индекс как "документ" — и агент начнёт отвечать на основе своих же старых реплик вместо содержимого настоящих файлов, что напрямую противоречит требованию про приоритет источника (Шаг 4).

`agent/memory.py`:

```python
import json, pathlib
import ollama

SESSIONS_DIR = pathlib.Path("sessions")
SESSIONS_DIR.mkdir(exist_ok=True)

SUMMARY_EVERY = 5   # суммаризировать каждые 5 пар "пользователь-агент"
KEEP_TAIL_PAIRS = 2 # столько последних пар оставляем дословно


def save_history(session_id: str, history: list):
    path = SESSIONS_DIR / f"{session_id}.jsonl"
    path.write_text("\n".join(json.dumps(m, ensure_ascii=False) for m in history))


def load_history(session_id: str) -> list:
    path = SESSIONS_DIR / f"{session_id}.jsonl"
    if not path.exists():
        return []
    return [json.loads(line) for line in path.read_text().splitlines() if line]


def _user_turns(history: list) -> int:
    return sum(1 for m in history if m["role"] == "user")


def maybe_summarize(history: list, model: str) -> list:
    """Раз в SUMMARY_EVERY пар реплик сжимает старую часть диалога в резюме,
    оставляя последние KEEP_TAIL_PAIRS пар дословно."""
    system = [m for m in history if m["role"] == "system"]
    rest = [m for m in history if m["role"] != "system"]

    if _user_turns(rest) < SUMMARY_EVERY:
        return history

    tail_len = KEEP_TAIL_PAIRS * 2  # пара = user + assistant (+ возможные tool-сообщения)
    to_summarize, tail = rest[:-tail_len], rest[-tail_len:]

    convo_text = "\n".join(f"{m['role']}: {m.get('content', '')}" for m in to_summarize)
    resp = ollama.chat(model=model, messages=[
        {"role": "system", "content": (
            "Сожми диалог ниже в резюме на 5-8 предложений: ключевые факты, "
            "принятые решения, договорённости. Ничего не добавляй от себя."
        )},
        {"role": "user", "content": convo_text},
    ])
    summary_msg = {
        "role": "system",
        "content": f"[Резюме предыдущей части диалога]: {resp['message']['content']}",
    }
    return system + [summary_msg] + tail
```

Суммаризация запускается автоматически в конце `run_agent` (Шаг 7) — пользователю не нужно ничего делать вручную. Это отдельный, "дешёвый" вызов модели, использующий тот же локальный Ollama, так что офлайн-режим не нарушается.

**Продолжение сессии.** По умолчанию (`agent "..."`) каждый запуск — новая сессия с одноразовым `SESSION_ID`, и после выхода она сохраняется в `sessions/<id>.jsonl`, но не подхватывается автоматически. Чтобы продолжить именно её или именованный проект:

```bash
agent --continue "..."                 # продолжить последнюю сессию (sessions/last.jsonl)
agent --session project-x "..."        # своя сессия под конкретную задачу/проект
```

## Шаг 9. Консольный вызов одной командой

```bash
cd ~/local-agent
source .venv/bin/activate   # если ещё не активировано с предыдущего шага
pip install -e .
```

`pyproject.toml` (минимум):
```toml
[project]
name = "local-agent"
version = "0.1.0"

[project.scripts]
agent = "agent.main:app"
```

`pip install -e .` кладёт исполняемый файл `agent` внутрь окружения — `.venv/bin/agent` (на Windows — `.venv\Scripts\agent.exe`). Отсюда два разных сценария использования:

**Вариант А — через активацию (по умолчанию).** Команда `agent` доступна, только пока `.venv` активировано в текущем терминале:
```bash
source ~/local-agent/.venv/bin/activate
agent "объясни, что в файле lecture3.pdf написано про производные"
agent            # интерактивный режим
```
В новом терминале, где `.venv` не активировано, `agent` не найдётся (`command not found`) — придётся активировать заново.

**Вариант Б — глобальный вызов без активации.** Чтобы `agent` работала сразу после открытия терминала, одной командой, без `source ...`, сделайте симлинк на исполняемый файл из окружения:
```bash
mkdir -p ~/.local/bin   # если ещё нет
ln -s ~/local-agent/.venv/bin/agent ~/.local/bin/agent
```
Убедитесь, что `~/.local/bin` есть в `PATH` (обычно уже добавлен в `~/.bashrc`/`~/.zshrc`; если нет — добавьте `export PATH="$HOME/.local/bin:$PATH"`). После этого:
```bash
agent "..."      # работает из любого нового терминала без активации venv
```
Symlink не копирует файл, а лишь указывает на него — python-интерпретатор и все пакеты (`ollama`, `chromadb` и т.д.) при запуске всё равно берутся из `.venv`, потому что путь до интерпретатора зашит в шебанг (`#!`) самого скрипта `agent`. Деактивировать `.venv` (`deactivate`) в этом варианте вообще не требуется — он не мешает работе симлинка.

## Шаг 10. Интеграция с VSCode

Для написания и правки кода прямо в IDE проще не встраивать это в собственный агент, а подключить готовый инструмент к тому же локальному Ollama-серверу:

- **Cline** (расширение VS Code) — открытый агент с правкой файлов и выполнением команд, поддерживает локальные модели через Ollama-эндпоинт `http://localhost:11434`.
- **Aider** — CLI-агент, который редактирует файлы и коммитит изменения в git; удобно запускать во встроенном терминале VSCode: `aider --model ollama/qwen2.5-coder:14b-instruct-q4_K_M`.
- Встроенный в VS Code **Language Model Provider** — с недавних версий можно подключить Ollama напрямую в чат Copilot без сторонних расширений (пункт в настройках "Manage Models" → Ollama). Проверьте актуальность в свежей документации VS Code, так как эта возможность появилась относительно недавно.

Ваш собственный CLI-агент (`agent`) при этом остаётся отдельным инструментом для скиллов/RAG/офлайн-браузера — VSCode-расширения такой логикой обычно не обладают.

## Проверка результата

- [ ] `ollama list` показывает скачанные модели без интернета.
- [ ] `agent "какой сейчас час"` вызывает скилл `get_time`.
- [ ] `agent "создай скилл, который переводит число из десятичной в двоичную систему"` — модель вызывает `create_skill`, в `skills/` появляется новая папка, и следующий вызов уже пользуется новым скиллом.
- [ ] После заливки двух документов в `docs/` и вызова `ingest()` для обоих: вопрос "по файлу А.pdf..." не подмешивает данные из Б.pdf (проверяется по `where`-фильтру в Chroma).
- [ ] При выключенной сети всё работает, кроме `browse_web`, который явно сообщает об ошибке сети.
- [ ] `agent "напиши и выполни код на C++, который..."` возвращает вывод из Docker-песочницы, а не из локального интерпретатора.
- [ ] После 5 пар реплик в `history` появляется системное сообщение `[Резюме предыдущей части диалога]`, а не полный текст старых сообщений.
- [ ] `agent --continue` после закрытия консоли подхватывает `sessions/last.jsonl` и помнит контекст предыдущего запуска.
- [ ] `agent "создай скилл, который умеет ходить в интернет по requests"` — агент не создаёт скилл молча и не отказывает без объяснений, а сообщает конкретную причину (`импорт модуля 'requests' запрещён...`) и ждёт подтверждения `[y/N]` в терминале.

## Что дорабатывать дальше

- Более умный роутинг "явно указанного документа" (сейчас — по точному совпадению имени файла; можно добавить fuzzy-matching или алиасы в конфиге).
- Логирование одобренных исключений `create_skill` (кто, когда и какое ограничение снял) — сейчас подтверждение через `input()` нигде не фиксируется, кроме самого факта появления небезопасного скилла на диске.
- Суммаризация сейчас использует ту же основную модель — при частых длинных диалогах может быть удобнее держать для этого отдельную маленькую модель (`qwen2.5:3b`), чтобы не грузить основную.
- Если 14B модели работают медленно на CPU — попробуйте `qwen2.5:7b-instruct-q4_K_M` или GGUF-квантование Q4_K_S.
