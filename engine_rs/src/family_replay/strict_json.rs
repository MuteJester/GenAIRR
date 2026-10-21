//! Duplicate-key rejecting bounded JSON parser for the new artifact only.
use super::FamilyReplayLimits;
use serde::de::{self, DeserializeSeed, MapAccess, SeqAccess, Visitor};
use serde_json::{Map, Number, Value};
use std::cell::Cell;
use std::fmt;

struct Seed<'a> {
    limits: &'a FamilyReplayLimits,
    count: &'a Cell<u64>,
    depth: u64,
}
impl<'de> DeserializeSeed<'de> for Seed<'_> {
    type Value = Value;
    fn deserialize<D: de::Deserializer<'de>>(self, d: D) -> Result<Value, D::Error> {
        if self.depth > self.limits.max_depth {
            return Err(de::Error::custom("resource limit: JSON depth"));
        }
        let n = self
            .count
            .get()
            .checked_add(1)
            .ok_or_else(|| de::Error::custom("JSON value count overflow"))?;
        if n > self.limits.max_values {
            return Err(de::Error::custom("resource limit: JSON values"));
        }
        self.count.set(n);
        d.deserialize_any(self)
    }
}
impl<'de> Visitor<'de> for Seed<'_> {
    type Value = Value;
    fn expecting(&self, f: &mut fmt::Formatter) -> fmt::Result {
        f.write_str("bounded JSON value")
    }
    fn visit_unit<E: de::Error>(self) -> Result<Value, E> {
        Ok(Value::Null)
    }
    fn visit_none<E: de::Error>(self) -> Result<Value, E> {
        Ok(Value::Null)
    }
    fn visit_bool<E: de::Error>(self, v: bool) -> Result<Value, E> {
        Ok(Value::Bool(v))
    }
    fn visit_i64<E: de::Error>(self, v: i64) -> Result<Value, E> {
        Ok(Value::Number(v.into()))
    }
    fn visit_u64<E: de::Error>(self, v: u64) -> Result<Value, E> {
        Ok(Value::Number(v.into()))
    }
    fn visit_f64<E: de::Error>(self, v: f64) -> Result<Value, E> {
        Number::from_f64(v)
            .map(Value::Number)
            .ok_or_else(|| E::custom("nonfinite JSON number"))
    }
    fn visit_str<E: de::Error>(self, v: &str) -> Result<Value, E> {
        if v.len() as u64 > self.limits.max_string_bytes {
            return Err(E::custom("resource limit: JSON string"));
        }
        Ok(Value::String(v.to_owned()))
    }
    fn visit_string<E: de::Error>(self, v: String) -> Result<Value, E> {
        if v.len() as u64 > self.limits.max_string_bytes {
            return Err(E::custom("resource limit: JSON string"));
        }
        Ok(Value::String(v))
    }
    fn visit_seq<A: SeqAccess<'de>>(self, mut seq: A) -> Result<Value, A::Error> {
        let mut values = Vec::new();
        while let Some(v) = seq.next_element_seed(Seed {
            limits: self.limits,
            count: self.count,
            depth: self.depth + 1,
        })? {
            values.push(v);
        }
        Ok(Value::Array(values))
    }
    fn visit_map<A: MapAccess<'de>>(self, mut map: A) -> Result<Value, A::Error> {
        let mut values = Map::new();
        while let Some(key) = map.next_key::<String>()? {
            if key.len() as u64 > self.limits.max_string_bytes {
                return Err(de::Error::custom("resource limit: JSON key"));
            }
            if values.contains_key(&key) {
                return Err(de::Error::custom(format!("duplicate JSON key: {key}")));
            }
            let value = map.next_value_seed(Seed {
                limits: self.limits,
                count: self.count,
                depth: self.depth + 1,
            })?;
            values.insert(key, value);
        }
        Ok(Value::Object(values))
    }
}
pub fn parse(text: &str, limits: &FamilyReplayLimits) -> Result<Value, String> {
    if text.len() as u64 > limits.max_file_bytes {
        return Err("resource limit: file bytes".into());
    }
    // serde_json's generic Value parser otherwise promotes an overflowing JSON
    // integer to f64. Preserve the declared integer/float distinction before
    // numeric conversion; quoted decimal/hex strings are unaffected.
    let bytes = text.as_bytes();
    let (mut i, mut quoted, mut escaped) = (0, false, false);
    while i < bytes.len() {
        let b = bytes[i];
        if quoted {
            if escaped {
                escaped = false;
            } else if b == b'\\' {
                escaped = true;
            } else if b == b'"' {
                quoted = false;
            }
            i += 1;
            continue;
        }
        if b == b'"' {
            quoted = true;
            i += 1;
            continue;
        }
        if b == b'-' || b.is_ascii_digit() {
            let start = i;
            i += 1;
            while i < bytes.len()
                && (bytes[i].is_ascii_digit()
                    || matches!(bytes[i], b'.' | b'e' | b'E' | b'+' | b'-'))
            {
                i += 1;
            }
            let token = &text[start..i];
            if !token.bytes().any(|b| matches!(b, b'.' | b'e' | b'E')) {
                let valid = if token.starts_with('-') {
                    token.parse::<i64>().is_ok()
                } else {
                    token.parse::<u64>().is_ok()
                };
                if !valid {
                    return Err("JSON integer outside [-2**63,2**64-1]".into());
                }
            }
            continue;
        }
        i += 1;
    }
    let count = Cell::new(0);
    let mut de = serde_json::Deserializer::from_str(text);
    let value = Seed {
        limits,
        count: &count,
        depth: 0,
    }
    .deserialize(&mut de)
    .map_err(|e| e.to_string())?;
    de.end().map_err(|e| e.to_string())?;
    Ok(value)
}
pub fn from_str<T: serde::de::DeserializeOwned>(
    text: &str,
    limits: &FamilyReplayLimits,
) -> Result<T, String> {
    serde_json::from_value(parse(text, limits)?).map_err(|e| e.to_string())
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn rejects_duplicate_nested_key_and_limits() {
        let l = FamilyReplayLimits::default();
        assert!(parse("{\"a\":{\"x\":1,\"x\":2}}", &l)
            .unwrap_err()
            .contains("duplicate"));
        assert!(parse(
            "[[[0]]]",
            &FamilyReplayLimits {
                max_depth: 2,
                ..l.clone()
            }
        )
        .is_err());
        assert!(parse(
            "[1,2,3]",
            &FamilyReplayLimits {
                max_values: 3,
                ..l.clone()
            }
        )
        .is_err());
        assert!(parse(
            "\"abcd\"",
            &FamilyReplayLimits {
                max_string_bytes: 3,
                ..l
            }
        )
        .is_err());
        assert!(parse("18446744073709551616", &FamilyReplayLimits::default()).is_err());
    }
}
