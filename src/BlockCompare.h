#if !defined(__BLOCKCOMPARATOR_H_)
#define __BLOCKCOMPARATOR_H_

template <typename T> class BlockComparator {
  const T *ptr;
  const size_t size;

public:
  BlockComparator(T *ptr, size_t size) : ptr(ptr), size(size) {}
  /**
   * Method to compare two blocks.
   * It relies on comparing features and statistical data
   * of each parameter. At the end, it returns a notion of
   * which block has more energy than the other.
   *
   * @param lhs: ponter to the left block.
   * @param rhs: pointer to the right block.
   * @param size: length of both arrays
   * @return
   *     -1 if lhs is smaller than rhs  <br/>
   *      0 if lhs equal rhs            <br/>
   *      1 if lhs is bigger than rhs   <br/>
   *
   */
  static int compare(const T *lhs, const T *rhs, const size_t size);

  /* Operator overloads */
  inline bool operator==(const T *other) {
    return compare(this->ptr, other) == 0;
  }
  inline bool operator!=(const T *other) {
    return compare(this->ptr, other) != 0;
  }
  inline bool operator<(const T *other) {
    return compare(this->ptr, other) < 0;
  }
  inline bool operator>(const T *other) {
    return compare(this->ptr, other) > 0;
  }
  inline bool operator<=(const T *other) {
    return compare(this->ptr, other) <= 0;
  }
  inline bool operator>=(const T *other) {
    return compare(this->ptr, other) >= 0;
  }
};

template <typename T>
int BlockComparator<T>::compare(const T *lhs, const T *rhs, const size_t size) {
  // Código aqui.
  // O tamanho dos arrays é dado pela variavel size.
  // Podes usar ela para fazer os loops, caso necessário.
  return 0;
}

#endif // __BLOCKCOMPARATOR_H_