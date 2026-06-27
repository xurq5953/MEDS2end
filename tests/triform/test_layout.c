#include "test_common.h"

static void test_layout_for_n(int n)
{
  Fq M[MEDS_n * MEDS_n * MEDS_n];
  Fq selector[MEDS_n];
  Fq selected[MEDS_n * MEDS_n];

  for (size_t i = 0; i < triform_element_count(n); i++)
    M[i] = (Fq)((i + 1u) % MEDS_p);

  CHECK(triform_slice_size(n) == (size_t)n * (size_t)n);
  CHECK(triform_element_count(n) == (size_t)n * (size_t)n * (size_t)n);

  for (int slice = 0; slice < n; slice++)
  {
    CHECK(triform_slice(M, slice, n) == M + (size_t)slice * triform_slice_size(n));
    CHECK(triform_slice_const(M, slice, n) == M + (size_t)slice * triform_slice_size(n));

    for (int row = 0; row < n; row++)
      for (int col = 0; col < n; col++)
      {
        size_t idx = triform_index(slice, row, col, n);
        CHECK(idx == (size_t)slice * triform_slice_size(n) + (size_t)row * (size_t)n + (size_t)col);
        CHECK(M[idx] == triform_slice_const(M, slice, n)[row * n + col]);
      }

    fill_zero(selector, (size_t)n);
    selector[slice] = 1;
    triform_matrix_at_w(selected, M, selector, n);
    expect_mat_equal(selected, triform_slice_const(M, slice, n), n);
  }
}

int main(void)
{
  test_layout_for_n(1);
  test_layout_for_n(2);
  test_layout_for_n(3);
  test_layout_for_n(MEDS_n);
  return 0;
}
