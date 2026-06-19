namespace ybi 
{

inline uint64_t MurmurHash64A(const uint8_t *key, size_t len, uint64_t seed)
{
    const uint64_t m = 0xc6a4a7935bd1e995ull;
    const int r = 47;

    uint64_t h = seed ^ (len * m);

    const uint8_t *end = key + 8 * (len / 8);

    while (key != end)
    {
        uint64_t k;
        std::memcpy(&k, key, sizeof(uint64_t));
        key += 8;

        k *= m;
        k ^= k >> r;
        k *= m;

        h ^= k;
        h *= m;
    }

    switch (len & 7)
    {
        case 7: h ^= uint64_t(key[6]) << 48;
        case 6: h ^= uint64_t(key[5]) << 40;
        case 5: h ^= uint64_t(key[4]) << 32;
        case 4: h ^= uint64_t(key[3]) << 24;
        case 3: h ^= uint64_t(key[2]) << 16;
        case 2: h ^= uint64_t(key[1]) << 8;
        case 1: h ^= uint64_t(key[0]); h *= m;
    };

    h ^= h >> r;
    h *= m;
    h ^= h >> r;

    return h;
}

}
