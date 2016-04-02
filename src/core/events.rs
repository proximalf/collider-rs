// Copyright 2016 Matthew D. Michelotti
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//   http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

use std::collections::BTreeMap;
use std::cmp::Ordering;
use std::f64;
use std::hash::{Hash, Hasher};
use core::HitboxId;
use util::{TightSet, n64, N64, OneOrTwo};

#[derive(Copy, Clone)]
pub struct EventKey {
    time: N64,
    index: u64
}

impl EventKey {
    fn time(&self) -> f64 {
        self.time.into()
    }
}

impl PartialEq for EventKey {
    fn eq(&self, rhs: &EventKey) -> bool {
        self.index == rhs.index
    }
}

impl Eq for EventKey { }

impl Hash for EventKey {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.index.hash(state)
    }
}

impl PartialOrd for EventKey {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for EventKey {
    fn cmp(&self, other: &Self) -> Ordering {
        if self.time == other.time {
            self.index.cmp(&other.index)
        } else {
            self.time.cmp(&other.time)
        }
    }
}

pub trait EventKeysMap {
    fn event_keys_mut(&mut self, id: HitboxId) -> &mut TightSet<EventKey>;
}

#[derive(Copy, Clone)]
pub enum InternalEvent {
    PanicSmallHitbox(HitboxId),
    PanicDurationPassed(HitboxId),
    Reiterate(HitboxId),
    Collide(HitboxId, HitboxId),
    Separate(HitboxId, HitboxId)
}

impl InternalEvent {
    fn other_id(self, id: HitboxId) -> Option<HitboxId> {
        self.involved_hitbox_ids().other_id(id)
    }
    
    fn involved_hitbox_ids(self) -> OneOrTwo<HitboxId> {
        match self {
            InternalEvent::PanicSmallHitbox(id) | InternalEvent::PanicDurationPassed(id)
                | InternalEvent::Reiterate(id) => OneOrTwo::One(id),
            InternalEvent::Collide(a, b) | InternalEvent::Separate(a, b) => OneOrTwo::Two(a, b)
        }
    }
}

pub struct EventManager {
    events: BTreeMap<EventKey, InternalEvent>,
    next_event_index: u64
}

impl EventManager {
    pub fn new() -> EventManager {
        EventManager { events : BTreeMap::new(), next_event_index : 0 }
    }

    pub fn add_solitaire_event(&mut self, time: f64, event: InternalEvent, key_set: &mut TightSet<EventKey>) {
        if let Some(key) = self.new_event_key(time) {
            assert!(self.events.insert(key, event).is_none(), "illegal state");
            assert!(key_set.insert(key), "illegal state");
        }
    }
    
    pub fn add_pair_event(&mut self, time: f64, event: InternalEvent, first_key_set: &mut TightSet<EventKey>,
        second_key_set: &mut TightSet<EventKey>)
    {
        if let Some(key) = self.new_event_key(time) {
            assert!(self.events.insert(key, event).is_none(), "illegal state");
            assert!(first_key_set.insert(key), "illegal state");
            assert!(second_key_set.insert(key), "illegal state");
        }
    }
    
    pub fn clear_related_events<M: EventKeysMap>(&mut self, id: HitboxId, key_set: &mut TightSet<EventKey>, map: &mut M) {
        for key in key_set.iter() {
            let event = self.events.remove(key).unwrap();
            if let Some(other_id) = event.other_id(id) {
                assert!(map.event_keys_mut(other_id).remove(key), "illegal state");
            }
        }
        key_set.clear();
    }
    
    fn new_event_key(&mut self, time: f64) -> Option<EventKey> {
        if time == f64::INFINITY {
            None
        } else {
            let result = EventKey { time: n64(time), index: self.next_event_index };
            self.next_event_index += 1;
            Some(result)
        }
    }
    
    pub fn peek_time(&self) -> f64 {
        self.peek_key().map(|key| key.time()).unwrap_or(f64::INFINITY)
    }
    
    pub fn next<M: EventKeysMap>(&mut self, time: f64, map: &mut M) -> Option<InternalEvent> {
        if let Some(key) = self.peek_key() {
            if key.time() == time {
                let event = self.events.remove(&key).unwrap();
                for id in event.involved_hitbox_ids().iter() {
                    assert!(map.event_keys_mut(id).remove(&key), "illegal state");
                }
                Some(event)
            } else {
                None
            }
        } else {
            None
        }
    }
    
    fn peek_key(&self) -> Option<EventKey> {
        self.events.keys().next().map(|&key| key)
    }
}