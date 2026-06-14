namespace ybi
{
inline bool IsPowerOfTwo(size_t value)
{
    return value != 0 && (value & (value - 1)) == 0;
}
} // namespace ybi
