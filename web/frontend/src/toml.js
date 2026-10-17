// Thin wrapper around smol-toml for the browser's raw TOML editor.
// Parse and stringify TOML that matches the strict control schema shape.
import { parse, stringify } from "smol-toml";

export function parseToml(text) {
  return parse(text);
}

export function dumpToml(obj) {
  return stringify(obj);
}
