pub struct MinstdRand {
    state: u64,
}

impl MinstdRand {
    const A: u64 = 48271;
    const M: u64 = 2147483647; // 2^31 - 1

    pub fn new(seed: u32) -> Self {
        Self { state: seed as u64 }
    }

    pub fn get_state(&self) -> u64 {
        self.state
    }

    /// Returns next integer in [1, M-1]
    fn next_int(&mut self) -> u64 {
        self.state = (self.state * Self::A) % Self::M;
        self.state
    }

    /// Returns next double in [0, 1) matching boost::uniform_01
    pub fn next_double(&mut self) -> f64 {
        let val = self.next_int();
        (val as f64 - 1.0) / (Self::M as f64 - 1.0)
    }
}

/// Mersenne Twister 19937 matching boost::mt19937.
/// Used as the main RNG for distance picking and coordinate generation.
pub struct Mt19937 {
    state: [u32; 624],
    index: usize,
}

impl Mt19937 {
    const N: usize = 624;
    const M: usize = 397;

    pub fn new(seed: u32) -> Self {
        let mut mt = Mt19937 {
            state: [0u32; 624],
            index: 624,
        };
        mt.state[0] = seed;
        for i in 1..Self::N {
            mt.state[i] = 1812433253u32
                .wrapping_mul(mt.state[i - 1] ^ (mt.state[i - 1] >> 30))
                .wrapping_add(i as u32);
        }
        mt
    }

    fn twist(&mut self) {
        for i in 0..Self::N {
            let y = (self.state[i] & 0x80000000) | (self.state[(i + 1) % Self::N] & 0x7fffffff);
            self.state[i] = self.state[(i + Self::M) % Self::N] ^ (y >> 1);
            if y & 1 != 0 {
                self.state[i] ^= 0x9908b0df;
            }
        }
        self.index = 0;
    }

    /// Generate next u32 value
    pub fn next_u32(&mut self) -> u32 {
        if self.index >= Self::N {
            self.twist();
        }
        let mut y = self.state[self.index];
        self.index += 1;

        y ^= y >> 11;
        y ^= (y << 7) & 0x9d2c5680;
        y ^= (y << 15) & 0xefc60000;
        y ^= y >> 18;
        y
    }

    /// Generate next double in [0, 1) matching boost::uniform_real<>(0, 1)
    /// Formula: mt_output / (mt_max + 1.0) = mt_output / 4294967296.0
    pub fn next_double(&mut self) -> f64 {
        self.next_u32() as f64 / 4294967296.0
    }
}
