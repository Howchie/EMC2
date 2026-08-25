#!/usr/bin/env python3
"""Generate a deterministic, static C++/R compatibility inventory.

This intentionally does not load R or import a package.  Facts that require
running a clean R session are retained as ``unresolved`` records instead of
being guessed or dropped.
"""

from __future__ import annotations

import argparse
import bisect
import hashlib
import json
import os
import re
import subprocess
import sys
from pathlib import Path
from typing import Any, Iterable

SCHEMA_VERSION = "0.1"
TARGET_FILES = (
    "R/RcppExports.R",
    "src/RcppExports.cpp",
    "src/col_registry.h",
    "src/particle_ll.cpp",
    "NAMESPACE",
)


def read_text(root: Path, rel: str) -> str:
    path = root / rel
    return path.read_text(encoding="utf-8", errors="replace") if path.exists() else ""


def line_number(text: str, offset: int) -> int:
    return bisect.bisect_right([0] + [m.start() for m in re.finditer("\\n", text)], offset)


def source(rel: str, text: str, offset: int) -> dict[str, Any]:
    return {"path": rel, "line": line_number(text, offset)}


def strip_r_comments(text: str) -> str:
    """Keep strings while removing # comments for lightweight R scanning."""
    out: list[str] = []
    i = 0
    in_string: str | None = None
    escaped = False
    while i < len(text):
        ch = text[i]
        if in_string:
            out.append(ch)
            if escaped:
                escaped = False
            elif ch == "\\":
                escaped = True
            elif ch == in_string:
                in_string = None
            i += 1
            continue
        if ch in "'\"":
            in_string = ch
            out.append(ch)
            i += 1
        elif ch == "#":
            j = text.find("\n", i)
            if j < 0:
                out.append(" " * (len(text) - i))
                break
            out.extend(" " * (j - i))
            out.append("\n")
            i = j + 1
        else:
            out.append(ch)
            i += 1
    return "".join(out)


def scan_balanced(text: str, start: int, opening: str = "(", closing: str = ")") -> int:
    """Return the matching delimiter offset, or -1 on an incomplete expression."""
    depth = 0
    in_string: str | None = None
    escaped = False
    i = start
    while i < len(text):
        ch = text[i]
        if in_string:
            if escaped:
                escaped = False
            elif ch == "\\":
                escaped = True
            elif ch == in_string:
                in_string = None
            i += 1
            continue
        if ch == "/" and i + 1 < len(text) and text[i + 1] == "/":
            nl = text.find("\n", i + 2)
            i = len(text) if nl < 0 else nl + 1
            continue
        if ch == "/" and i + 1 < len(text) and text[i + 1] == "*":
            end_comment = text.find("*/", i + 2)
            i = len(text) if end_comment < 0 else end_comment + 2
            continue
        if ch in "'\"":
            in_string = ch
        elif ch == opening:
            depth += 1
        elif ch == closing:
            depth -= 1
            if depth == 0:
                return i
        i += 1
    return -1


def split_top_level(text: str, separator: str = ",") -> list[str]:
    parts: list[str] = []
    start = 0
    stack: list[str] = []
    in_string: str | None = None
    escaped = False
    pairs = {"(": ")", "[": "]", "{": "}"}
    for i, ch in enumerate(text):
        if in_string:
            if escaped:
                escaped = False
            elif ch == "\\":
                escaped = True
            elif ch == in_string:
                in_string = None
            continue
        if ch in "'\"":
            in_string = ch
        elif ch in pairs:
            stack.append(pairs[ch])
        elif stack and ch == stack[-1]:
            stack.pop()
        elif ch == separator and not stack:
            parts.append(text[start:i])
            start = i + 1
    parts.append(text[start:])
    return parts


def quoted_strings(text: str) -> list[str]:
    return [m.group(2) for m in re.finditer(r"(['\"])((?:\\.|(?!\1).)*)\1", text, re.S)]


def unescape_r_string(value: str) -> str:
    # The source facts used in names are ordinary ASCII literals.  Decode only
    # escapes needed for a stable display; avoid Python's unicode escape quirks.
    return re.sub(r"\\([\\'\"nrt])", lambda m: {"n": "\n", "r": "\r", "t": "\t"}.get(m.group(1), m.group(1)), value)


def function_records(text: str) -> list[dict[str, Any]]:
    clean = strip_r_comments(text)
    records: list[dict[str, Any]] = []
    pattern = re.compile(r"(?m)^[ \t]*([.]?[A-Za-z][A-Za-z0-9_.]*)[ \t]*(?:<-|=)[ \t]*function[ \t]*\(")
    for m in pattern.finditer(clean):
        name = m.group(1)
        open_paren = clean.find("(", m.start())
        close_paren = scan_balanced(clean, open_paren)
        if close_paren < 0:
            continue
        j = close_paren + 1
        while j < len(clean) and clean[j].isspace():
            j += 1
        if j >= len(clean) or clean[j] != "{":
            # One-expression functions are retained because exported model
            # constructors may delegate directly to a helper.
            end = clean.find("\n", j)
            if end < 0:
                end = len(text)
            records.append({
                "name": name,
                "start": m.start(),
                "params_start": open_paren + 1,
                "params_end": close_paren,
                "body_start": j,
                "body_end": end,
                "body": text[j:end],
                "header": text[m.start():close_paren + 1],
            })
            continue
        close_brace = scan_balanced(clean, j, "{", "}")
        if close_brace < 0:
            close_brace = len(text) - 1
        records.append({
            "name": name,
            "start": m.start(),
            "params_start": open_paren + 1,
            "params_end": close_paren,
            "body_start": j + 1,
            "body_end": close_brace,
            "body": text[j + 1:close_brace],
            "header": text[m.start():close_paren + 1],
        })
    return records


def parse_arguments(header: str) -> list[dict[str, str]]:
    open_paren = header.find("(")
    close_paren = header.rfind(")")
    if open_paren < 0 or close_paren < open_paren:
        return []
    out: list[dict[str, str]] = []
    for raw in split_top_level(header[open_paren + 1:close_paren]):
        item = raw.strip()
        if not item:
            continue
        bits = split_top_level(item, "=")
        name = bits[0].strip()
        default = "" if len(bits) == 1 else "=".join(bits[1:]).strip()
        out.append({"name": name, "default": default})
    return out


def expression_after(body: str, field: str) -> tuple[str, int] | None:
    m = re.search(r"(?m)^[ \t]*(?:[A-Za-z][A-Za-z0-9_.]*[ \t]*<-?[ \t]*)?" + re.escape(field) + r"[ \t]*(?:<-|=)[ \t]*", body)
    if not m:
        # c_name commonly appears as a list member with a leading comma.
        m = re.search(r"\b" + re.escape(field) + r"[ \t]*(?:<-|=)[ \t]*", body)
    if not m:
        return None
    start = m.end()
    i = start
    stack: list[str] = []
    in_string: str | None = None
    escaped = False
    while i < len(body):
        ch = body[i]
        if in_string:
            if escaped:
                escaped = False
            elif ch == "\\":
                escaped = True
            elif ch == in_string:
                in_string = None
        elif ch in "'\"":
            in_string = ch
        elif ch in "([{":
            stack.append({"(": ")", "[": "]", "{": "}"}[ch])
        elif stack and ch == stack[-1]:
            stack.pop()
        elif not stack and ch in ",\n":
            break
        i += 1
    return body[start:i].strip(), m.start()


def c_name_info(body: str) -> dict[str, Any] | None:
    found = expression_after(body, "c_name")
    if not found:
        return None
    expr, _ = found
    # Resolve the common local-variable form: c_name <- paste0(...); list(c_name = c_name).
    if re.fullmatch(r"c_name", expr):
        local = re.search(r"(?m)^[ \t]*c_name[ \t]*(?:<-|=)[ \t]*(.+)$", body)
        if local:
            expr = local.group(1).strip().rstrip(",")
    strings = sorted({unescape_r_string(s) for s in quoted_strings(expr) if s})
    static_values: list[str] = []
    if len(strings) == 1 and not re.search(r"\b(paste0?|ifelse|if|match\.arg)\b", expr):
        static_values = strings
    suffix_tokens = sorted({s for s in strings if s.startswith("_")})
    return {
        "expression": expr,
        "static_values": static_values,
        "tokens": strings,
        "suffix_tokens": suffix_tokens,
        "resolved_statically": bool(static_values),
    }


def vector_expression(body: str, field: str) -> str | None:
    found = expression_after(body, field)
    if not found:
        return None
    expr, _ = found
    return expr


def named_vector_names(expr: str | None) -> list[str]:
    if not expr:
        return []
    out: list[str] = []
    for m in re.finditer(r"(?:^|,|\\n)\s*(?:[\"']([^\"']+)[\"']|([A-Za-z][A-Za-z0-9_.]*))\s*=", expr):
        out.append(m.group(1) or m.group(2))
    return sorted(dict.fromkeys(out))


def plain_vector_names(expr: str | None) -> list[str]:
    if not expr:
        return []
    return sorted(dict.fromkeys(unescape_r_string(s) for s in quoted_strings(expr)))


def model_inventory(root: Path, unresolved: list[dict[str, Any]]) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[str]]:
    constructors: list[dict[str, Any]] = []
    model_files: list[dict[str, Any]] = []
    for path in sorted((root / "R").glob("model_*.R"), key=lambda p: p.as_posix()):
        rel = path.relative_to(root).as_posix()
        text = path.read_text(encoding="utf-8", errors="replace")
        funcs = function_records(text)
        found = []
        for rec in funcs:
            body = rec["body"]
            cinfo = c_name_info(body)
            # Dot-prefixed helpers can construct a model list, but are not
            # user-facing model constructors themselves.  Their expressions
            # remain available through the constructor records below.
            if not cinfo or rec["name"].startswith("."):
                continue
            p_expr = vector_expression(body, "p_types")
            canon_expr = vector_expression(body, "p_types_canonical")
            p_names = named_vector_names(p_expr)
            if not p_names:
                # A p_types=c("a",...) expression is uncommon but still useful.
                p_names = plain_vector_names(p_expr)
            canon_names = plain_vector_names(canon_expr)
            optional = sorted(set(p_names) - set(canon_names)) if p_names and canon_names else []
            item: dict[str, Any] = {
                "name": rec["name"],
                "kind": "constructor",
                "source": source(rel, text, rec["start"]),
                "arguments": parse_arguments(rec["header"]),
                "c_name": cinfo,
                "type_tokens": sorted(set(quoted_strings(vector_expression(body, "type") or ""))),
                "p_types": {
                    "expression": p_expr,
                    "names": p_names,
                    "canonical_expression": canon_expr,
                    "canonical_names": canon_names,
                    "optional_or_nuisance_names": optional,
                },
            }
            constructors.append(item)
            found.append(rec["name"])
            unresolved.append({
                "kind": "r_evaluation",
                "subject": f"model:{rec['name']}",
                "fact": "constructor metadata and generated c_name family",
                "reason": "requires evaluating the constructor in a clean R session; static expressions/tokens are included above",
                "source": source(rel, text, rec["start"]),
            })
        # Exported constructors that delegate to a helper (for example
        # BTAwL -> .btawl_constructor) have no literal c_name in their body.
        # Keep the constructor and mark the delegated family unresolved.
        for rec in funcs:
            if rec["name"].startswith(".") or rec["name"] in found:
                continue
            prior = text[max(0, rec["start"] - 2000):rec["start"]]
            body = rec["body"]
            exported = "@export" in prior
            model_like = bool(re.search(r"(?:constructor|\bp_types\b|\brfun\b|\bpfun\b|\blog_likelihood\b|\bc_name\b|\btype\s*=)", body))
            is_constructor_name = rec["name"][0].isupper() or rec["name"] in {"softmax", "hUVSD"}
            if not exported or not model_like or not is_constructor_name:
                continue
            item = {
                "name": rec["name"],
                "kind": "constructor",
                "source": source(rel, text, rec["start"]),
                "arguments": parse_arguments(rec["header"]),
                "c_name": {
                    "expression": "delegated constructor call; inspect helper at runtime",
                    "static_values": [],
                    "tokens": sorted(set(quoted_strings(body))),
                    "suffix_tokens": sorted({s for s in quoted_strings(body) if s.startswith("_")}),
                    "resolved_statically": False,
                },
                "p_types": {
                    "expression": None,
                    "names": [],
                    "canonical_expression": None,
                    "canonical_names": [],
                    "optional_or_nuisance_names": [],
                },
            }
            constructors.append(item)
            found.append(rec["name"])
            unresolved.append({
                "kind": "r_evaluation",
                "subject": f"model:{rec['name']}",
                "fact": "delegated constructor c_name and parameter metadata",
                "reason": "constructor delegates to a helper or branch that requires clean-R evaluation",
                "source": source(rel, text, rec["start"]),
            })
        # Include exported aliases such as GOMPERTZ <- GOM and Gompertz <- GOM.
        aliases = dict(re.findall(r"(?m)^\s*([A-Za-z][A-Za-z0-9_.]*)\s*<-\s*([A-Za-z][A-Za-z0-9_.]*)\s*$", text))
        known_names = {x["name"] for x in constructors}
        for alias, target in sorted(aliases.items()):
            if target in known_names and alias not in known_names:
                target_item = next(x for x in constructors if x["name"] == target)
                constructors.append({
                    "name": alias,
                    "kind": "alias",
                    "target": target,
                    "source": source(rel, text, text.find(alias)),
                    "c_name": target_item["c_name"],
                    "p_types": target_item["p_types"],
                })
                found.append(alias)
        has_model_metadata = bool(re.search(r"\bc_name\s*=", text))
        model_files.append({
            "path": rel,
            "represented": bool(found),
            "constructor_names": sorted(found),
            "contains_model_metadata": has_model_metadata,
        })
    constructors.sort(key=lambda x: (x["name"], x["source"]["path"], x["source"]["line"]))
    model_files.sort(key=lambda x: x["path"])
    return constructors, model_files, [x["name"] for x in constructors]


def rcpp_inventory(root: Path) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    rel = "R/RcppExports.R"
    text = read_text(root, rel)
    wrappers: list[dict[str, Any]] = []
    for rec in function_records(text):
        body = rec["body"]
        calls = list(re.finditer(r"\.Call\s*\(\s*[`\"]([^`\"]+)[`\"]", body))
        if not calls:
            continue
        for call in calls:
            symbol = call.group(1)
            wrappers.append({
                "r_name": rec["name"],
                "symbol": symbol,
                "arguments": parse_arguments(rec["header"]),
                "source": source(rel, text, rec["start"]),
                "call_source": source(rel, text, rec["body_start"] + call.start()),
            })
    cpp_rel = "src/RcppExports.cpp"
    cpp = read_text(root, cpp_rel)
    exports: list[dict[str, Any]] = []
    for m in re.finditer(r"(?m)^\s*RcppExport\s+SEXP\s+(_[A-Za-z0-9_]+)\s*\(([^)]*)\)", cpp):
        symbol = m.group(1)
        exports.append({"symbol": symbol, "arity": len([x for x in split_top_level(m.group(2)) if x.strip()]), "source": source(cpp_rel, cpp, m.start())})
    registrations: list[dict[str, Any]] = []
    for m in re.finditer(r'\{"(_[A-Za-z0-9_]+)",\s*\(DL_FUNC\)\s*&(_[A-Za-z0-9_]+),\s*(\d+)\}', cpp):
        registrations.append({"symbol": m.group(1), "target": m.group(2), "arity": int(m.group(3)), "source": source(cpp_rel, cpp, m.start())})
    wrappers.sort(key=lambda x: (x["r_name"], x["symbol"]))
    exports.sort(key=lambda x: x["symbol"])
    registrations.sort(key=lambda x: x["symbol"])
    reg_symbols = {x["symbol"] for x in registrations}
    exported_symbols = {x["symbol"] for x in exports}
    pairs = []
    for wrapper in wrappers:
        symbol = wrapper["symbol"]
        generated = next((x for x in exports if x["symbol"] == symbol), None)
        registered = next((x for x in registrations if x["symbol"] == symbol), None)
        pairs.append({
            "r_name": wrapper["r_name"],
            "symbol": symbol,
            "generated_wrapper": generated is not None,
            "registered": symbol in reg_symbols,
            "generated": generated,
            "registration": registered,
            "source": wrapper["source"],
        })
    return {
        "r_wrappers": wrappers,
        "generated_wrappers": exports,
        "registration": registrations,
        "pairs": pairs,
        "registration_enabled": bool(re.search(r"R_registerRoutines\s*\(", cpp) and re.search(r"R_useDynamicSymbols\s*\([^,]+,\s*FALSE", cpp)),
    }, []


def col_registry_inventory(root: Path) -> list[dict[str, Any]]:
    rel = "src/col_registry.h"
    text = read_text(root, rel)
    clean = re.sub(r"//.*", "", text)
    out: list[dict[str, Any]] = []
    parent_comment = ""
    for m in re.finditer(r"(?m)^\s*namespace\s+([A-Za-z_][A-Za-z0-9_]*)\s*\{", clean):
        ns = m.group(1)
        if ns == "emc2col":
            continue
        close = scan_balanced(clean, clean.find("{", m.start()), "{", "}")
        if close < 0:
            continue
        body = clean[m.end():close]
        enum_m = re.search(r"enum\s*:\s*int\s*\{(.*?)\}\s*;", body, re.S)
        enum_required: list[str] = []
        enum_optional: list[str] = []
        if enum_m:
            after = False
            for bit in split_top_level(enum_m.group(1)):
                token = bit.strip()
                if not token:
                    continue
                name = token.split("=", 1)[0].strip()
                if name == "N_REQ":
                    after = True
                elif after:
                    enum_optional.append(name)
                else:
                    enum_required.append(name)
        for spec_m in re.finditer(r"inline\s+ColSpec\s+(\w+)\s*\(\)\s*\{(.*?)\n\s*\}", body, re.S):
            spec_name = spec_m.group(1)
            sb = spec_m.group(2)
            array_m = re.search(r"static\s+const\s+char\*\s+n\[\]\s*=\s*\{(.*?)\};", sb, re.S)
            if not array_m:
                continue
            names = [unescape_r_string(x) for x in quoted_strings(array_m.group(1))]
            label_m = re.search(r"return\s*\{\s*n\s*,\s*N_REQ\s*,\s*[\"']([^\"']+)", sb)
            label = label_m.group(1) if label_m else ns
            # The nearest model comment is intentionally recovered from the
            # original text, not from preprocessed C++.
            before = text[:m.start()]
            cm = list(re.finditer(r"R/model_[A-Za-z0-9_]+\.R", before))
            model_file = cm[-1].group(0) if cm else None
            out.append({
                "namespace": ns,
                "spec": spec_name,
                "label": label,
                "required": names,
                "required_count": len(names),
                "optional": sorted(dict.fromkeys(enum_optional)),
                "enum_required": enum_required,
                "source": source(rel, text, m.start()),
                "model_file": model_file,
            })
    return sorted(out, key=lambda x: (x["namespace"], x["spec"]))


def dispatch_inventory(root: Path) -> dict[str, Any]:
    rel = "src/particle_ll.cpp"
    text = read_text(root, rel)
    marker = text.find("resolve_race_model_adapter")
    if marker < 0:
        return {"race_precedence": [], "route_symbols": [], "suffix_tokens": [], "conditions": []}
    body_start = text.find("{", marker)
    body_end = scan_balanced(text, body_start, "{", "}")
    body = text[body_start + 1:body_end if body_end >= 0 else len(text)]
    # Top-level branch starts in this resolver.  Nested checks use four or more
    # spaces and are deliberately excluded from precedence ordering.
    starts: list[tuple[int, str, str]] = []
    for m in re.finditer(r"(?m)^( {2})(?:(?:if)|(?:\}\s*else\s+if))\s*\((.*)$", body):
        line = m.group(2)
        if "type_std.find" not in line:
            continue
        cond = line
        # Include a continuation line for `GOM || type_std.find("GOMP")`.
        p = m.end()
        while "{" not in cond and p < len(body):
            nxt = body.find("\n", p)
            if nxt < 0:
                break
            line2_end = body.find("\n", nxt + 1)
            if line2_end < 0:
                line2_end = len(body)
            continuation = body[nxt + 1:line2_end].strip()
            cond += " " + continuation
            p = line2_end
            if "{" in continuation:
                break
        tokens = [unescape_r_string(x) for x in quoted_strings(cond) if "type_std.find" in cond]
        # Restrict tokens to arguments of type_std.find, avoiding constants in
        # branch comments/defaults.
        tokens = [unescape_r_string(x) for x in re.findall(r"type_std\.find\s*\(\s*[\"']([^\"']+)", cond)]
        starts.append((m.start(), cond.strip(), tokens[0] if tokens else ""))
    # A leading `if` is usually captured; de-duplicate continuation artifacts.
    unique: list[tuple[int, str, str]] = []
    seen = set()
    for item in starts:
        key = (item[0], item[2])
        if key not in seen:
            seen.add(key)
            unique.append(item)
    routes: list[dict[str, Any]] = []
    precedence: list[dict[str, Any]] = []
    conditions: set[str] = set()
    for i, (offset, cond, token) in enumerate(unique):
        end = unique[i + 1][0] if i + 1 < len(unique) else len(body)
        chunk = body[offset:end]
        tokens = sorted(set(re.findall(r"type_std\.find\s*\(\s*[\"']([^\"']+)", chunk[: min(len(chunk), chunk.find("{") + 1 if "{" in chunk else len(chunk))])))
        if not tokens:
            tokens = [token] if token else []
        conditions.update(tokens)
        if token == "GNG":
            # GNG is a feature flag on the DDM path, not a separate race
            # adapter branch.  Keep it in conditions but not route_symbols.
            conditions.update(tokens)
            continue
        symbols = sorted(set(re.findall(r"=\s*&([A-Za-z_][A-Za-z0-9_]*)", chunk)))
        specs = sorted(set(re.findall(r"emc2col::([A-Za-z_][A-Za-z0-9_]*)::(?:spec(?:_\w+)?)", chunk)))
        route = {"order": i, "condition": cond, "match_tokens": tokens, "symbols": symbols, "column_namespaces": specs, "source": source(rel, text, marker + body_start - marker + offset)}
        routes.append(route)
        if i and token != "GNG":
            precedence.append({"before": unique[i - 1][2], "after": token, "order": i, "reason": "ordered substring dispatch branch"})
    all_tokens = sorted(set(re.findall(r"type_std\.find\s*\(\s*[\"']([^\"']+)", body)))
    collisions = []
    branch_tokens = [r["match_tokens"][0] for r in routes if r["match_tokens"]]
    for earlier_i, earlier in enumerate(branch_tokens):
        for later in branch_tokens[earlier_i + 1:]:
            if earlier and later and (later in earlier or earlier in later):
                collisions.append({"before": earlier, "after": later, "reason": "substring collision requires specialized-before-generic review"})
    # Capture non-race dispatcher names and route symbols from the same source.
    nonrace = sorted(set(re.findall(r"(?:type_std\s*==\s*|type\s*==\s*)[\"']([^\"']+)[\"']", text)))
    return {
        "race_precedence": routes,
        "precedence_constraints": sorted(precedence + collisions, key=lambda x: (x.get("order", 10**9), x["before"], x["after"])),
        "route_symbols": sorted(set(s for r in routes for s in r["symbols"])),
        "column_namespaces": sorted(set(s for r in routes for s in r["column_namespaces"])),
        "suffix_tokens": all_tokens,
        "conditions": sorted(conditions),
        "nonrace_exact_types": nonrace,
        "source": source(rel, text, marker),
    }


def target_references(root: Path, unresolved: list[dict[str, Any]]) -> dict[str, Any]:
    refs: dict[str, list[dict[str, Any]]] = {"direct_call": [], "xptr": [], "custom_trend": [], "custom_kernel": [], "source_cpp": []}
    # Direct calls are kept precise for all R files; the requested generated R
    # file is still the principal source of wrapper calls.
    scan_paths: list[Path] = []
    for directory in (root / "R", root / "src"):
        if directory.exists():
            scan_paths.extend(p for p in directory.rglob("*") if p.is_file() and p.suffix in {".R", ".r", ".cpp", ".h", ".cc", ".hpp"})
    for path in sorted(scan_paths, key=lambda p: p.relative_to(root).as_posix()):
        rel = path.relative_to(root).as_posix()
        text = path.read_text(encoding="utf-8", errors="replace")
        for m in re.finditer(r"\.Call\s*\(\s*[`\"]([^`\"]+)[`\"]", text):
            refs["direct_call"].append({"symbol": m.group(1), "source": source(rel, text, m.start())})
        for m in re.finditer(r"(?:Rcpp::)?XPtr(?:\s*<\s*([^>]+)>\s*)?", text):
            refs["xptr"].append({"type": (m.group(1) or "unspecified").strip(), "source": source(rel, text, m.start())})
        for pattern, key in ((r"custom[_ ]trend|funptr", "custom_trend"), (r"custom[_ ]kernel|register_kernel", "custom_kernel"), (r"sourceCpp", "source_cpp")):
            for m in re.finditer(pattern, text, re.I):
                refs[key].append({"match": m.group(0), "context": text[m.start():text.find("\n", m.start()) if text.find("\n", m.start()) >= 0 else len(text)].strip(), "source": source(rel, text, m.start())})
    for key in refs:
        unique = {(json.dumps(x, sort_keys=True)): x for x in refs[key]}
        refs[key] = sorted(unique.values(), key=lambda x: (x["source"]["path"], x["source"]["line"], json.dumps(x, sort_keys=True)))
    if not refs["xptr"]:
        unresolved.append({"kind": "ffi_scan", "subject": "Rcpp::XPtr", "fact": "external-pointer ownership and finalizers", "reason": "no XPtr expression was found in the scanned source files"})
    if not refs["source_cpp"]:
        unresolved.append({"kind": "ffi_scan", "subject": "Rcpp::sourceCpp", "fact": "runtime compilation users", "reason": "no sourceCpp call was found in scanned source files"})
    return refs


def namespace_inventory(root: Path) -> dict[str, Any]:
    rel = "NAMESPACE"
    text = read_text(root, rel)
    exports = sorted(set(re.findall(r"(?m)^\s*export\(([^)]+)\)", text)))
    imports = sorted(set(re.findall(r"(?m)^\s*importFrom\(([^,]+),\s*([^)]+)\)", text)))
    dyn = sorted(set(re.findall(r"(?m)^\s*useDynLib\(([^)]+)\)", text)))
    return {"exports": exports, "imports": [{"package": a.strip(), "name": b.strip()} for a, b in imports], "use_dynlib": dyn, "source": {"path": rel}}


def revision(root: Path) -> str:
    # Revision is provenance, not an input to the manifest.  Unknown is valid
    # in source archives and keeps generated output reproducible there.
    try:
        proc = subprocess.run(["git", "-C", str(root), "rev-parse", "HEAD"], capture_output=True, text=True, check=True, timeout=2)
        value = proc.stdout.strip()
        return value if re.fullmatch(r"[0-9a-fA-F]{7,64}", value) else "unknown"
    except (OSError, subprocess.SubprocessError):
        return "unknown"


def make_manifest(root: Path) -> tuple[dict[str, Any], list[str]]:
    unresolved: list[dict[str, Any]] = []
    constructors, model_files, _ = model_inventory(root, unresolved)
    rcpp, _ = rcpp_inventory(root)
    registry = col_registry_inventory(root)
    dispatch = dispatch_inventory(root)
    references = target_references(root, unresolved)
    namespace = namespace_inventory(root)
    files: list[dict[str, Any]] = []
    for rel in TARGET_FILES:
        path = root / rel
        if path.exists():
            data = path.read_bytes()
            files.append({"path": rel, "sha256": hashlib.sha256(data).hexdigest(), "bytes": len(data)})
        else:
            unresolved.append({"kind": "source_file", "subject": rel, "fact": "required inventory input", "reason": "file is absent from this source tree"})
    # Explicitly record clean-R facts not derivable from syntax.
    unresolved.extend([
        {"kind": "provenance", "subject": "build_profile", "fact": "compiler/flags/OpenMP profile", "reason": "build profile is not represented in the requested source inputs"},
        {"kind": "r_evaluation", "subject": "all_model_constructors", "fact": "legal argument combinations, transforms, bounds, and runtime routes", "reason": "requires evaluating R constructors in a clean R session"},
    ])
    checks: dict[str, Any] = {"errors": [], "warnings": []}
    r_symbols = {x["symbol"] for x in rcpp["r_wrappers"]}
    registration_symbols = {x["symbol"] for x in rcpp["registration"]}
    missing_registration = sorted(r_symbols - registration_symbols)
    if missing_registration:
        checks["errors"].append({"kind": "missing_registration", "symbols": missing_registration, "message": "RcppExports.R .Call wrapper(s) are absent from src/RcppExports.cpp registration"})
    unrepresented = [x["path"] for x in model_files if x["contains_model_metadata"] and not x["represented"]]
    if unrepresented:
        checks["errors"].append({"kind": "unrepresented_model_files", "files": unrepresented, "message": "model_*.R file contains c_name metadata but no discovered constructor"})
    if not rcpp["registration_enabled"]:
        checks["warnings"].append({"kind": "registration_mode", "message": "Rcpp registration could not be proven from generated C++"})
    # Pair/registration mismatch is actionable but generated standalone helpers
    # are expected and therefore remain warnings, not check failures.
    extra_registration = sorted(registration_symbols - r_symbols)
    if extra_registration:
        checks["warnings"].append({"kind": "standalone_generated_entries", "symbols": extra_registration, "message": "generated registration entries have no RcppExports.R wrapper; they may be direct C++ helpers"})
    unresolved.sort(key=lambda x: (x.get("kind", ""), x.get("subject", ""), json.dumps(x, sort_keys=True)))
    manifest: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "inventory": "EMC2 C++ compatibility inventory",
        "provenance": {
            "revision": revision(root),
            "build_profile": "unknown",
            "execution": "static standard-library-only scan; no R evaluation",
            "files": files,
        },
        "models": {"constructors": constructors, "files": model_files},
        "rcpp": rcpp,
        "columns": {"registry": registry},
        "dispatch": dispatch,
        "references": references,
        "namespace": namespace,
        "unresolved": unresolved,
        "checks": checks,
    }
    return manifest, [json.dumps(x, sort_keys=True) for x in checks["errors"]]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd(), help="package root (default: current directory)")
    parser.add_argument("--output", "-o", type=Path, help="write JSON to this path (default: stdout)")
    parser.add_argument("--check", action="store_true", help="validate Rcpp wrapper registration and model-file coverage")
    args = parser.parse_args(argv)
    root = args.root.resolve()
    manifest, errors = make_manifest(root)
    payload = json.dumps(manifest, indent=2, sort_keys=True, ensure_ascii=False) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload, encoding="utf-8")
    else:
        sys.stdout.write(payload)
    if args.check and errors:
        sys.stderr.write("compatibility inventory check failed:\n")
        for error in manifest["checks"]["errors"]:
            sys.stderr.write("- " + error["message"] + ": " + ", ".join(error.get("symbols", error.get("files", []))) + "\n")
        return 1
    if args.check:
        sys.stderr.write("compatibility inventory check passed (no missing Rcpp registrations or model-file coverage gaps)\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
