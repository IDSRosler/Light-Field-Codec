# `Block.h` API

**Block** is a template wrapper to any contiguous memory location. It allows to index either single or multidimension arrays.
Because **Block** is a wrapper, no ownership is taken regarding the underlying memory space.

## Usage

#### Create a block

To construct `Block` objects, a pointer and a lenght has to be passed to the constructor.

```cpp
int sarray[120];            // Stack allocated array
int *harray = new int[120]; // Heap allocated array

Block bs(sarray, 120);
Block bh(harray, 120);
```

#### Reshaping

You may reshape the indexer to any valid shape. A shape is valid if the product of the dimensions remains the same.

```cpp
bs.reshape(5, 24);      // valid shape: 5 * 24 == 120
bs.reshape(5, 4, 6);    // valid shape: 5 * 4 * 6 == 120
bs.reshape(3, 30);      // invalid shape: 3 * 30 == 90 != 120
                        // throws std::length_error
```

You may call `shape()` to get the current block shape. It returns a `std::vector<std::size_t>` object.

```cpp
bs.reshape(5, 24);
auto shape = bs.shape();
// shape[0] == 5
// shape[1] == 24
// shape.size() == 2
```

#### Index

To index the block, simply pass the values to the operator `()`

```cpp
bs.reshape(5, 24);      // 2D array
int x = bs(3, 10);

bs.reshape(5, 4, 6);    // 3D array
int y = bs(3, 2, 5);
```

It is also possible to assign values.

```cpp
bs.reshape(5, 24);
bs(1, 3) = 45;          // Valid
```

Please note that going beyond the length of a dimension is not allowed.

```cpp
bs.reshape(120)
int invalid_value = bs(120);    // throws std::out_of_range
bs.reshape(5, 24);
int other_invalid = bs(5, 10);  // throws std::out_of_range

```

Because `reshape` allows to add _fake_ dimensons, indexing should change as well.

```cpp
bs.reshape(1, 120, 1);
int value = bs(0, 12, 0);
```

#### Views

A view is a subset of a block. Because the unerlying memory locations are never copied, all assignment done on a view will be an assignment done on the block itself.
To create views, use the `Range` operator. It takes two parameters `start` and `end`.

```cpp
using Range = Block<int>::Range;
bs.reshape(5, 24);
auto v = bs.view({Range(0, 3), Range(5, 20)});  // Note the curly braces.
// v(1, 1) == bs(1, 6)
```

It is allowed to call `shape()` on a view.

```cpp
...
auto v = bs.view({Range(0, 3), Range(5, 20)});  // Note the curly braces.
auto shape = v.shape();
// shape = {3, 15}
```

#### Iterators

There are 4 iterators options available.

##### 1) `flat_view()` and `flat_size()`

```cpp
auto v = bs.flat_view();
for (int i = 0; i < v.flat_size(); i++) {
    // use v(i) to reference the object
}
```

##### 2) Forward Iterator `begin()` and `end()`

```cpp
for(auto it = bs.begin(); it != bs.end(); it++) {
    // use *it to reference the object
}
```

##### 3) Reverse Iterator `rbegin()` and `rend()`

```cpp
for(auto it = bs.rbegin(); it != bs.rend(); it++) {
    // use *it to reference the object
}
```

##### 4) For range syntax

```cpp
for (auto &it: bs) {
    // use it to reference the object
}
```
