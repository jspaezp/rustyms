//! Benchmark-only Rust allocator counters. No allocation or locking in callbacks.
use std::{
    alloc::{GlobalAlloc, Layout, System},
    sync::atomic::{AtomicUsize, Ordering::Relaxed},
};

pub(super) struct CountingAllocator {
    allocs: AtomicUsize,
    frees: AtomicUsize,
    reallocs: AtomicUsize,
    requested: AtomicUsize,
    growth: AtomicUsize,
    realloc_bytes: AtomicUsize,
    live: AtomicUsize,
    peak: AtomicUsize,
}
#[derive(Debug)]
#[allow(dead_code)] // Fields are emitted through Debug outside the measured region.
pub(super) struct Report {
    pub(super) alloc_calls: usize,
    pub(super) free_calls: usize,
    pub(super) realloc_calls: usize,
    pub(super) alloc_requested_bytes: usize,
    pub(super) realloc_requested_bytes: usize,
    pub(super) growth_bytes: usize,
    pub(super) baseline_live_bytes: usize,
    pub(super) end_live_bytes: usize,
    pub(super) peak_live_bytes: usize,
}
pub(super) struct Start(usize, usize, usize, usize, usize, usize, usize);
impl CountingAllocator {
    pub(super) const fn new() -> Self {
        Self {
            allocs: AtomicUsize::new(0),
            frees: AtomicUsize::new(0),
            reallocs: AtomicUsize::new(0),
            requested: AtomicUsize::new(0),
            growth: AtomicUsize::new(0),
            realloc_bytes: AtomicUsize::new(0),
            live: AtomicUsize::new(0),
            peak: AtomicUsize::new(0),
        }
    }
    // requested includes fresh allocations and full realloc destinations;
    // growth includes fresh allocations and only positive realloc size deltas.
    fn grow(&self, bytes: usize) {
        self.growth.fetch_add(bytes, Relaxed);
        let live = self.live.fetch_add(bytes, Relaxed) + bytes;
        self.peak.fetch_max(live, Relaxed);
    }
    pub(super) fn begin(&self) -> Start {
        let live = self.live.load(Relaxed);
        self.peak.store(live, Relaxed);
        Start(
            self.allocs.load(Relaxed),
            self.frees.load(Relaxed),
            self.reallocs.load(Relaxed),
            self.requested.load(Relaxed),
            self.growth.load(Relaxed),
            live,
            self.realloc_bytes.load(Relaxed),
        )
    }
    pub(super) fn finish(&self, start: Start) -> Report {
        let realloc_bytes = self.realloc_bytes.load(Relaxed) - start.6;
        Report {
            alloc_calls: self.allocs.load(Relaxed) - start.0,
            free_calls: self.frees.load(Relaxed) - start.1,
            realloc_calls: self.reallocs.load(Relaxed) - start.2,
            alloc_requested_bytes: self.requested.load(Relaxed) - start.3 - realloc_bytes,
            realloc_requested_bytes: realloc_bytes,
            growth_bytes: self.growth.load(Relaxed) - start.4,
            baseline_live_bytes: start.5,
            end_live_bytes: self.live.load(Relaxed),
            peak_live_bytes: self.peak.load(Relaxed),
        }
    }
}
// SAFETY: all operations forward the original pointer/layout to System. Counters
// do not allocate, dereference pointers, or change allocation layout/ownership.
unsafe impl GlobalAlloc for CountingAllocator {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        let ptr = unsafe { System.alloc(layout) };
        if !ptr.is_null() {
            self.allocs.fetch_add(1, Relaxed);
            self.requested.fetch_add(layout.size(), Relaxed);
            self.grow(layout.size());
        }
        ptr
    }
    unsafe fn alloc_zeroed(&self, layout: Layout) -> *mut u8 {
        let ptr = unsafe { System.alloc_zeroed(layout) };
        if !ptr.is_null() {
            self.allocs.fetch_add(1, Relaxed);
            self.requested.fetch_add(layout.size(), Relaxed);
            self.grow(layout.size());
        }
        ptr
    }
    unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
        self.frees.fetch_add(1, Relaxed);
        self.live.fetch_sub(layout.size(), Relaxed);
        unsafe { System.dealloc(ptr, layout) };
    }
    unsafe fn realloc(&self, ptr: *mut u8, layout: Layout, size: usize) -> *mut u8 {
        let result = unsafe { System.realloc(ptr, layout, size) };
        if !result.is_null() {
            self.reallocs.fetch_add(1, Relaxed);
            self.requested.fetch_add(size, Relaxed);
            self.realloc_bytes.fetch_add(size, Relaxed);
            if size >= layout.size() {
                self.grow(size - layout.size());
            } else {
                self.live.fetch_sub(layout.size() - size, Relaxed);
            }
        }
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn tracks_growth_shrink_zeroed_and_free() {
        // An independent instance avoids test-harness/global allocator traffic.
        let allocator = CountingAllocator::new();
        let start = allocator.begin();
        unsafe {
            let small = Layout::from_size_align(16, 8).unwrap();
            let large = Layout::from_size_align(64, 8).unwrap();
            let ptr = allocator.alloc(small);
            assert!(!ptr.is_null());
            let ptr = allocator.realloc(ptr, small, 64);
            assert!(!ptr.is_null());
            let ptr = allocator.realloc(ptr, large, 16);
            assert!(!ptr.is_null());
            allocator.dealloc(ptr, small);
            let ptr = allocator.alloc_zeroed(small);
            assert!(!ptr.is_null());
            assert_eq!(std::slice::from_raw_parts(ptr, 16), &[0; 16]);
            allocator.dealloc(ptr, small);
        }
        let report = allocator.finish(start);
        assert_eq!(report.alloc_calls, 2);
        assert_eq!(report.free_calls, 2);
        assert_eq!(report.realloc_calls, 2);
        assert_eq!(report.alloc_requested_bytes, 32);
        assert_eq!(report.realloc_requested_bytes, 80);
        assert_eq!(report.growth_bytes, 80);
        assert_eq!(report.baseline_live_bytes, 0);
        assert_eq!(report.end_live_bytes, 0);
        assert_eq!(report.peak_live_bytes, 64);
    }
}
