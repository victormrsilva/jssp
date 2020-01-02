\ Model chvatal
\ LP format - for model browsing. Use MPS format to capture full model detail.
Maximize
  - 0.0001 u(ord(0,1)) - 0.0001 u(ord(0,2)) - 0.0001 u(ord(1,2))
   - 0.0001 u(ord(2,1)) - 0.0001 u(ord(2,2)) - 0.0001 u(triangle2(0,1,0))
   - 0.0001 u(triangle2(0,1,1)) - 0.0001 u(triangle2(0,1,2))
   - 0.0001 u(triangle2(0,2,0)) - 0.0001 u(triangle2(0,2,1))
   - 0.0001 u(triangle2(0,2,2)) - 0.0001 u(triangle2(1,2,0))
   - 0.0001 u(triangle2(1,2,1)) - 0.0001 u(triangle2(1,2,2))
   - 0.0001 u(MCLinV(0,1,0)) - 0.0001 u(MCLinV(1,0,0))
   - 0.0001 u(MCLin1(0,1,0)) - 0.0001 u(MCLin2(0,1,0))
   - 0.0001 u(MCLin4(0,1,0)) - 0.0001 u(MCLin8(1,0,0))
   - 0.0001 u(MCLin10(1,0,0)) - 0.0001 u(MCLinV(0,1,1))
   - 0.0001 u(MCLinV(1,0,1)) - 0.0001 u(MCLin2(0,1,1))
   - 0.0001 u(MCLin7(1,0,1)) - 0.0001 u(MCLin8(1,0,1))
   - 0.0001 u(MCLinV(0,1,2)) - 0.0001 u(MCLinV(1,0,2))
   - 0.0001 u(MCLin2(0,1,2)) - 0.0001 u(MCLin8(1,0,2))
   - 0.0001 u(MCLinV(0,2,0)) - 0.0001 u(MCLinV(2,0,0))
   - 0.0001 u(MCLin3(0,2,0)) - 0.0001 u(MCLin5(0,2,0))
   - 0.0001 u(MCLin10(2,0,0)) - 0.0001 u(MCLinV(0,2,1))
   - 0.0001 u(MCLinV(2,0,1)) - 0.0001 u(MCLin2(0,2,1))
   - 0.0001 u(MCLin3(0,2,1)) - 0.0001 u(MCLin7(2,0,1))
   - 0.0001 u(MCLinV(0,2,2)) - 0.0001 u(MCLinV(2,0,2))
   - 0.0001 u(MCLin3(0,2,2)) - 0.0001 u(MCLin5(0,2,2))
   - 0.0001 u(MCLin10(2,0,2)) - 0.0001 u(MCLinV(1,2,0))
   - 0.0001 u(MCLinV(2,1,0)) - 0.0001 u(MCLin2(1,2,0))
   - 0.0001 u(MCLin7(2,1,0)) - 0.0001 u(MCLinV(1,2,1))
   - 0.0001 u(MCLinV(2,1,1)) - 0.0001 u(MCLin2(1,2,1))
   - 0.0001 u(MCLin10(2,1,1)) - 0.0001 u(MCLinV(1,2,2))
   - 0.0001 u(MCLinV(2,1,2)) - 0.0001 u(MCLin1(1,2,2))
   - 0.0001 u(MCLin2(1,2,2)) - 0.0001 u(MCLin4(1,2,2))
   - 0.0001 u(MCLin8(2,1,2)) - 0.0001 u(MCLin10(2,1,2))
   - 0.0001 u(makespan(0)) - 0.0001 u(makespan(1)) - 0.0001 u(makespan(2))
   - 0.0001 u(cut_late_job_best(012,0,1,2))
   - 0.0001 u(cut_basic_best(012,1)) - 0.0001 u(cut_half_best(1,012,1))
   - 0.0001 u(cut_late_job_best(012,1,1,0))
   - 0.0001 u(cut_late_job_best(012,1,1,1))
   - 0.0001 u(cut_late_job_best(012,1,1,2))
   - 0.0001 u(cut_basic_best(012,2)) - 0.0001 u(cut_half_best(0,012,2))
   - 0.0001 u(cut_late_job_best(012,2,0,0))
   - 0.0001 u(cut_late_job_best(012,2,0,1))
   - 0.0001 u(cut_late_job_best(012,2,0,2)) + obj
Subject To
 objective: 21.3 a(C) + 6.3 a(x(0,1)) + 13.3 a(x(0,2)) + 20.3 a(x(0,0))
   + 11.3 a(x(1,0)) + 4.2 a(x(1,1)) + 14.3 a(x(1,2)) + 6.3 a(x(2,1))
   + 13.3 a(x(2,2)) + 19.3 a(x(2,0)) + 0.168 a(y(0,1,1))
   + 0.832 a(y(1,0,1)) + 0.1 a(y(0,1,2)) + 0.9 a(y(1,0,2))
   + 0.7136363636363637 a(y(0,2,1)) + 0.2863636363636363 a(y(2,0,1))
   + 0.3363636363636366 a(y(1,2,0)) + 0.6636363636363634 a(y(2,1,0))
   + 0.568 a(y(1,2,1)) + 0.432 a(y(2,1,1)) + 21.3 a(v(0,1,0))
   + 14.3 a(v(1,0,0)) + 14.3 a(u(1,0,0)) + 13.3 a(v(0,1,1))
   + 8.2 a(v(1,0,1)) + 1.176 a(u(0,1,1)) + 3.328 a(u(1,0,1))
   + 20.3 a(v(0,1,2)) + 21.3 a(v(1,0,2)) + 1.4 a(u(0,1,2))
   + 17.7 a(u(1,0,2)) + 21.3 a(v(0,2,0)) + 21.3 a(v(2,0,0))
   + 21.3 a(u(0,2,0)) + 6.3 a(u(2,0,0)) + 13.3 a(v(0,2,1))
   + 13.3 a(v(2,0,1)) + 4.99545454545455 a(u(0,2,1))
   + 2.00454545454545 a(u(2,0,1)) + 20.3 a(v(0,2,2)) + 19.3 a(v(2,0,2))
   + 20.3 a(u(0,2,2)) + 6.3 a(u(2,0,2)) + 14.3 a(v(1,2,0))
   + 21.3 a(v(2,1,0)) + 2.354545454545473 a(u(1,2,0))
   + 9.954545454545439 a(u(2,1,0)) + 8.2 a(v(1,2,1)) + 13.3 a(v(2,1,1))
   + 2.272 a(u(1,2,1)) + 9.324 a(u(2,1,1)) + 21.3 a(v(1,2,2))
   + 19.3 a(v(2,1,2)) + 19.3 a(u(2,1,2)) - a0 - obj = -0
 f_C: - u(makespan(0)) - u(makespan(1)) - u(makespan(2)) - a(C) - f(C)
   = -0
 f_x(0,1): u(ord(0,1)) + u(MCLinV(0,1,1)) + u(MCLinV(0,2,1))
   - 0.25 u(cut_clique(012,1)) - 7 u(cut_basic_best(012,1))
   - 7 u(cut_two_jobs(0,2,1)) - u(cut_half_best(0,012,1))
   - u(cut_late_job_best(012,1,0,0)) - u(cut_late_job_best(012,1,0,1))
   - u(cut_late_job_best(012,1,0,2))
   - 7 u(cut_basic_plus_epsilon_best(02,1,1))
   - 7 u(cut_basic_plus_epsilon_best(01,1,2)) - a(x(0,1)) - f(x(0,1)) = -0
 f_x(0,2): - u(ord(0,1)) + u(ord(0,2)) + u(MCLinV(0,1,2))
   + u(MCLinV(0,2,2)) - 7 u(cut_basic_best(012,2))
   - 7 u(cut_two_jobs(0,1,2)) - u(cut_half_best(0,012,2))
   - u(cut_late_job_best(012,2,0,0)) - u(cut_late_job_best(012,2,0,1))
   - u(cut_late_job_best(012,2,0,2))
   - 7 u(cut_basic_plus_epsilon_best(02,2,1))
   - 7 u(cut_basic_plus_epsilon_best(01,2,2)) - a(x(0,2)) - f(x(0,2)) = -0
 f_x(0,0): - u(ord(0,2)) + u(MCLinV(0,1,0)) + u(MCLinV(0,2,0))
   + u(makespan(0)) - u(cut_basic_best(02,0)) - 2 u(cut_two_jobs(0,2,0))
   - u(cut_half_best(0,02,0)) - u(cut_late_job_best(02,0,0,0))
   - u(cut_late_job_best(012,0,0,2))
   - u(cut_basic_plus_epsilon_best(01,0,2)) - a(x(0,0)) - f(x(0,0)) = -0
 f_x(1,0): - u(ord(1,1)) + u(ord(1,2)) + u(MCLinV(1,0,0))
   + u(MCLinV(1,2,0)) - u(cut_half_best(1,12,0))
   - u(cut_late_job_best(012,0,1,0)) - u(cut_late_job_best(012,0,1,1))
   - u(cut_late_job_best(012,0,1,2))
   - 3 u(cut_basic_plus_epsilon_best(01,0,2)) - a(x(1,0)) - f(x(1,0)) = -0
 f_x(1,1): u(ord(1,1)) + u(MCLinV(1,0,1)) + u(MCLinV(1,2,1))
   - 0.14286 u(cut_clique(012,1)) - 4 u(cut_basic_best(012,1))
   - 4 u(cut_basic_plus_epsilon_best(12,1,0)) - u(cut_half_best(1,012,1))
   - u(cut_late_job_best(012,1,1,0)) - u(cut_late_job_best(012,1,1,1))
   - u(cut_late_job_best(012,1,1,2))
   - 4 u(cut_basic_plus_epsilon_best(01,1,2)) - a(x(1,1)) - f(x(1,1)) = -0
 f_x(1,2): - u(ord(1,2)) + u(MCLinV(1,0,2)) + u(MCLinV(1,2,2))
   + u(makespan(1)) - 7 u(cut_basic_best(012,2)) - 7 u(cut_two_jobs(0,1,2))
   - 7 u(cut_basic_plus_epsilon_best(12,2,0)) - u(cut_half_best(1,012,2))
   - u(cut_late_job_best(012,2,1,0)) - u(cut_late_job_best(012,2,1,1))
   - u(cut_late_job_best(012,2,1,2))
   - 7 u(cut_basic_plus_epsilon_best(01,2,2)) - a(x(1,2)) - f(x(1,2)) = -0
 f_x(2,1): u(ord(2,1)) + u(MCLinV(2,0,1)) + u(MCLinV(2,1,1))
   - 7 u(cut_basic_best(012,1)) - 7 u(cut_two_jobs(0,2,1))
   - 7 u(cut_basic_plus_epsilon_best(12,1,0))
   - 7 u(cut_basic_plus_epsilon_best(02,1,1)) - a(x(2,1)) - f(x(2,1)) = -0
 f_x(2,2): - u(ord(2,1)) + u(ord(2,2)) + u(MCLinV(2,0,2))
   + u(MCLinV(2,1,2)) - 6 u(cut_basic_best(012,2))
   - 6 u(cut_basic_plus_epsilon_best(12,2,0))
   - 6 u(cut_basic_plus_epsilon_best(02,2,1)) - a(x(2,2)) - f(x(2,2)) = -0
 f_x(2,0): - u(ord(2,2)) + u(MCLinV(2,0,0)) + u(MCLinV(2,1,0))
   + u(makespan(2)) - 2 u(cut_basic_best(02,0)) - u(cut_two_jobs(0,2,0))
   - a(x(2,0)) - f(x(2,0)) = -0
 f_y(0,1,1): u(triangle2(0,1,1)) + 7 u(MCLin2(0,1,1)) + 30 u(MCLin3(0,1,1))
   - 30 u(MCLin4(0,1,1)) - 7 u(MCLin5(0,1,1)) - 29 u(MCLin9(1,0,1))
   + 4 u(MCLin10(1,0,1)) + 7 u(cut_half_best(1,012,1))
   + 7 u(cut_late_job_best(012,1,1,0)) + 7 u(cut_late_job_best(012,1,1,1))
   + 7 u(cut_late_job_best(012,1,1,2)) - a(y(0,1,1)) - f(y(0,1,1)) = -0
 f_y(1,0,1): u(triangle2(0,1,1)) + 4 u(MCLin7(1,0,1)) + 29 u(MCLin8(1,0,1))
   + 4 u(cut_half_best(0,012,1)) + 4 u(cut_late_job_best(012,1,0,0))
   + 4 u(cut_late_job_best(012,1,0,1)) + 4 u(cut_late_job_best(012,1,0,2))
   - a(y(1,0,1)) - f(y(1,0,1)) = -0
 f_y(0,1,2): u(triangle2(0,1,2)) + 14 u(MCLin2(0,1,2))
   + 37 u(MCLin3(0,1,2)) - 37 u(MCLin4(0,1,2)) - 14 u(MCLin5(0,1,2))
   - 36 u(MCLin9(1,0,2)) + 14 u(MCLin10(1,0,2))
   + 7 u(cut_half_best(1,012,2)) + 7 u(cut_late_job_best(012,2,1,0))
   + 7 u(cut_late_job_best(012,2,1,1)) + 7 u(cut_late_job_best(012,2,1,2))
   - a(y(0,1,2)) - f(y(0,1,2)) = -0
 f_y(1,0,2): u(triangle2(0,1,2)) + 14 u(MCLin7(1,0,2))
   + 36 u(MCLin8(1,0,2)) + 7 u(cut_half_best(0,012,2))
   + 7 u(cut_late_job_best(012,2,0,0)) + 7 u(cut_late_job_best(012,2,0,1))
   + 7 u(cut_late_job_best(012,2,0,2)) - a(y(1,0,2)) - f(y(1,0,2)) = -0
 f_y(0,2,1): u(triangle2(0,2,1)) + 7 u(MCLin2(0,2,1)) + 29 u(MCLin3(0,2,1))
   - 29 u(MCLin4(0,2,1)) - 7 u(MCLin5(0,2,1)) - 29 u(MCLin9(2,0,1))
   + 7 u(MCLin10(2,0,1)) - a(y(0,2,1)) - f(y(0,2,1)) = -0
 f_y(2,0,1): u(triangle2(0,2,1)) + 7 u(MCLin7(2,0,1)) + 29 u(MCLin8(2,0,1))
   + 7 u(cut_half_best(0,012,1)) + 7 u(cut_late_job_best(012,1,0,0))
   + 7 u(cut_late_job_best(012,1,0,1)) + 7 u(cut_late_job_best(012,1,0,2))
   - a(y(2,0,1)) - f(y(2,0,1)) = -0
 f_y(1,2,0): u(triangle2(1,2,0)) + 7 u(MCLin2(1,2,0)) + 37 u(MCLin3(1,2,0))
   - 37 u(MCLin4(1,2,0)) - 7 u(MCLin5(1,2,0)) - 34 u(MCLin9(2,1,0))
   + 15 u(MCLin10(2,1,0)) - 9 u(cut_late_job_best(012,0,0,2))
   - 9 u(cut_late_job_best(012,0,1,2))
   - 36 u(cut_basic_plus_epsilon_best(01,0,2)) - a(y(1,2,0)) - f(y(1,2,0))
   = -0
 f_y(2,1,0): u(triangle2(1,2,0)) + 15 u(MCLin7(2,1,0))
   + 34 u(MCLin8(2,1,0)) + 2 u(cut_half_best(1,12,0))
   + 2 u(cut_late_job_best(012,0,1,0)) + 2 u(cut_late_job_best(012,0,1,1))
   + 2 u(cut_late_job_best(012,0,1,2)) - a(y(2,1,0)) - f(y(2,1,0)) = -0
 f_y(1,2,1): u(triangle2(1,2,1)) + 4 u(MCLin2(1,2,1)) + 29 u(MCLin3(1,2,1))
   - 29 u(MCLin4(1,2,1)) - 4 u(MCLin5(1,2,1)) - 30 u(MCLin9(2,1,1))
   + 7 u(MCLin10(2,1,1)) - a(y(1,2,1)) - f(y(1,2,1)) = -0
 f_y(2,1,1): u(triangle2(1,2,1)) + 7 u(MCLin7(2,1,1)) + 30 u(MCLin8(2,1,1))
   + 7 u(cut_half_best(1,012,1)) + 7 u(cut_late_job_best(012,1,1,0))
   + 7 u(cut_late_job_best(012,1,1,1)) + 7 u(cut_late_job_best(012,1,1,2))
   - a(y(2,1,1)) - f(y(2,1,1)) = -0
 f_v(0,1,0): - u(MCLinV(0,1,0)) + u(MCLin3(0,1,0)) - u(MCLin5(0,1,0))
   - a(v(0,1,0)) - f(v(0,1,0)) = -0
 f_v(1,0,0): - u(MCLinV(1,0,0)) + u(MCLin8(1,0,0)) - u(MCLin10(1,0,0))
   - a(v(1,0,0)) - f(v(1,0,0)) = -0
 f_u(1,0,0): - u(MCLin6(1,0,0)) - u(MCLin7(1,0,0)) - u(MCLin8(1,0,0))
   + u(MCLin9(1,0,0)) + u(MCLin10(1,0,0)) - a(u(1,0,0)) - f(u(1,0,0)) = -0
 f_v(0,1,1): - u(MCLinV(0,1,1)) + u(MCLin3(0,1,1)) - u(MCLin5(0,1,1))
   - a(v(0,1,1)) - f(v(0,1,1)) = -0
 f_v(1,0,1): - u(MCLinV(1,0,1)) + u(MCLin8(1,0,1)) - u(MCLin10(1,0,1))
   - a(v(1,0,1)) - f(v(1,0,1)) = -0
 f_u(0,1,1): - u(MCLin1(0,1,1)) - u(MCLin2(0,1,1)) - u(MCLin3(0,1,1))
   + u(MCLin4(0,1,1)) + u(MCLin5(0,1,1)) - a(u(0,1,1)) - f(u(0,1,1)) = -0
 f_u(1,0,1): - u(MCLin6(1,0,1)) - u(MCLin7(1,0,1)) - u(MCLin8(1,0,1))
   + u(MCLin9(1,0,1)) + u(MCLin10(1,0,1)) - a(u(1,0,1)) - f(u(1,0,1)) = -0
 f_v(0,1,2): - u(MCLinV(0,1,2)) + u(MCLin3(0,1,2)) - u(MCLin5(0,1,2))
   - a(v(0,1,2)) - f(v(0,1,2)) = -0
 f_v(1,0,2): - u(MCLinV(1,0,2)) + u(MCLin8(1,0,2)) - u(MCLin10(1,0,2))
   - a(v(1,0,2)) - f(v(1,0,2)) = -0
 f_u(0,1,2): - u(MCLin1(0,1,2)) - u(MCLin2(0,1,2)) - u(MCLin3(0,1,2))
   + u(MCLin4(0,1,2)) + u(MCLin5(0,1,2)) - a(u(0,1,2)) - f(u(0,1,2)) = -0
 f_u(1,0,2): - u(MCLin6(1,0,2)) - u(MCLin7(1,0,2)) - u(MCLin8(1,0,2))
   + u(MCLin9(1,0,2)) + u(MCLin10(1,0,2)) - a(u(1,0,2)) - f(u(1,0,2)) = -0
 f_v(0,2,0): - u(MCLinV(0,2,0)) + u(MCLin3(0,2,0)) - u(MCLin5(0,2,0))
   - a(v(0,2,0)) - f(v(0,2,0)) = -0
 f_v(2,0,0): - u(MCLinV(2,0,0)) + u(MCLin8(2,0,0)) - u(MCLin10(2,0,0))
   - a(v(2,0,0)) - f(v(2,0,0)) = -0
 f_u(0,2,0): - u(MCLin1(0,2,0)) - u(MCLin2(0,2,0)) - u(MCLin3(0,2,0))
   + u(MCLin4(0,2,0)) + u(MCLin5(0,2,0)) - a(u(0,2,0)) - f(u(0,2,0)) = -0
 f_u(2,0,0): - u(MCLin6(2,0,0)) - u(MCLin7(2,0,0)) - u(MCLin8(2,0,0))
   + u(MCLin9(2,0,0)) + u(MCLin10(2,0,0)) - a(u(2,0,0)) - f(u(2,0,0)) = -0
 f_v(0,2,1): - u(MCLinV(0,2,1)) + u(MCLin3(0,2,1)) - u(MCLin5(0,2,1))
   - a(v(0,2,1)) - f(v(0,2,1)) = -0
 f_v(2,0,1): - u(MCLinV(2,0,1)) + u(MCLin8(2,0,1)) - u(MCLin10(2,0,1))
   - a(v(2,0,1)) - f(v(2,0,1)) = -0
 f_u(0,2,1): - u(MCLin1(0,2,1)) - u(MCLin2(0,2,1)) - u(MCLin3(0,2,1))
   + u(MCLin4(0,2,1)) + u(MCLin5(0,2,1)) - a(u(0,2,1)) - f(u(0,2,1)) = -0
 f_u(2,0,1): - u(MCLin6(2,0,1)) - u(MCLin7(2,0,1)) - u(MCLin8(2,0,1))
   + u(MCLin9(2,0,1)) + u(MCLin10(2,0,1)) - a(u(2,0,1)) - f(u(2,0,1)) = -0
 f_v(0,2,2): - u(MCLinV(0,2,2)) + u(MCLin3(0,2,2)) - u(MCLin5(0,2,2))
   - a(v(0,2,2)) - f(v(0,2,2)) = -0
 f_v(2,0,2): - u(MCLinV(2,0,2)) + u(MCLin8(2,0,2)) - u(MCLin10(2,0,2))
   - a(v(2,0,2)) - f(v(2,0,2)) = -0
 f_u(0,2,2): - u(MCLin1(0,2,2)) - u(MCLin2(0,2,2)) - u(MCLin3(0,2,2))
   + u(MCLin4(0,2,2)) + u(MCLin5(0,2,2)) - a(u(0,2,2)) - f(u(0,2,2)) = -0
 f_u(2,0,2): - u(MCLin6(2,0,2)) - u(MCLin7(2,0,2)) - u(MCLin8(2,0,2))
   + u(MCLin9(2,0,2)) + u(MCLin10(2,0,2)) - a(u(2,0,2)) - f(u(2,0,2)) = -0
 f_v(1,2,0): - u(MCLinV(1,2,0)) + u(MCLin3(1,2,0)) - u(MCLin5(1,2,0))
   - a(v(1,2,0)) - f(v(1,2,0)) = -0
 f_v(2,1,0): - u(MCLinV(2,1,0)) + u(MCLin8(2,1,0)) - u(MCLin10(2,1,0))
   - a(v(2,1,0)) - f(v(2,1,0)) = -0
 f_u(1,2,0): - u(MCLin1(1,2,0)) - u(MCLin2(1,2,0)) - u(MCLin3(1,2,0))
   + u(MCLin4(1,2,0)) + u(MCLin5(1,2,0)) - a(u(1,2,0)) - f(u(1,2,0)) = -0
 f_u(2,1,0): - u(MCLin6(2,1,0)) - u(MCLin7(2,1,0)) - u(MCLin8(2,1,0))
   + u(MCLin9(2,1,0)) + u(MCLin10(2,1,0)) - a(u(2,1,0)) - f(u(2,1,0)) = -0
 f_v(1,2,1): - u(MCLinV(1,2,1)) + u(MCLin3(1,2,1)) - u(MCLin5(1,2,1))
   - a(v(1,2,1)) - f(v(1,2,1)) = -0
 f_v(2,1,1): - u(MCLinV(2,1,1)) + u(MCLin8(2,1,1)) - u(MCLin10(2,1,1))
   - a(v(2,1,1)) - f(v(2,1,1)) = -0
 f_u(1,2,1): - u(MCLin1(1,2,1)) - u(MCLin2(1,2,1)) - u(MCLin3(1,2,1))
   + u(MCLin4(1,2,1)) + u(MCLin5(1,2,1)) - a(u(1,2,1)) - f(u(1,2,1)) = -0
 f_u(2,1,1): - u(MCLin6(2,1,1)) - u(MCLin7(2,1,1)) - u(MCLin8(2,1,1))
   + u(MCLin9(2,1,1)) + u(MCLin10(2,1,1)) - a(u(2,1,1)) - f(u(2,1,1)) = -0
 f_v(1,2,2): - u(MCLinV(1,2,2)) + u(MCLin3(1,2,2)) - u(MCLin5(1,2,2))
   - a(v(1,2,2)) - f(v(1,2,2)) = -0
 f_v(2,1,2): - u(MCLinV(2,1,2)) + u(MCLin8(2,1,2)) - u(MCLin10(2,1,2))
   - a(v(2,1,2)) - f(v(2,1,2)) = -0
 f_u(2,1,2): - u(MCLin6(2,1,2)) - u(MCLin7(2,1,2)) - u(MCLin8(2,1,2))
   + u(MCLin9(2,1,2)) + u(MCLin10(2,1,2)) - a(u(2,1,2)) - f(u(2,1,2)) = -0
 f: - 7 u(ord(0,1)) - 7 u(ord(0,2)) - 4 u(ord(1,1)) - 3 u(ord(1,2))
   - 7 u(ord(2,1)) - 6 u(ord(2,2)) + u(triangle2(0,1,1))
   + u(triangle2(0,1,2)) + u(triangle2(0,2,1)) + u(triangle2(1,2,0))
   + u(triangle2(1,2,1)) - u(MCLinV(0,1,0)) - 3 u(MCLinV(1,0,0))
   + 34 u(MCLin3(0,1,0)) - 15 u(MCLin5(0,1,0)) - 7 u(MCLin7(1,0,0))
   + 37 u(MCLin9(1,0,0)) - 7 u(MCLinV(0,1,1)) - 4 u(MCLinV(1,0,1))
   + 30 u(MCLin3(0,1,1)) - 7 u(MCLin5(0,1,1)) + 29 u(MCLin8(1,0,1))
   + 29 u(MCLin9(1,0,1)) - 7 u(MCLinV(0,1,2)) - 7 u(MCLinV(1,0,2))
   + 37 u(MCLin3(0,1,2)) - 14 u(MCLin5(0,1,2)) + 36 u(MCLin8(1,0,2))
   + 36 u(MCLin9(1,0,2)) - u(MCLinV(0,2,0)) - 2 u(MCLinV(2,0,0))
   - 15 u(MCLin2(0,2,0)) + 42 u(MCLin4(0,2,0)) + 43 u(MCLin8(2,0,0))
   + 86 u(MCLin9(2,0,0)) - 15 u(MCLin10(2,0,0)) - 7 u(MCLinV(0,2,1))
   - 7 u(MCLinV(2,0,1)) + 29 u(MCLin3(0,2,1)) - 7 u(MCLin5(0,2,1))
   + 29 u(MCLin8(2,0,1)) + 29 u(MCLin9(2,0,1)) - 7 u(MCLinV(0,2,2))
   - 6 u(MCLinV(2,0,2)) - 14 u(MCLin2(0,2,2)) + 36 u(MCLin4(0,2,2))
   + 36 u(MCLin8(2,0,2)) + 72 u(MCLin9(2,0,2)) - 13 u(MCLin10(2,0,2))
   - 3 u(MCLinV(1,2,0)) - 2 u(MCLinV(2,1,0)) + 37 u(MCLin3(1,2,0))
   - 7 u(MCLin5(1,2,0)) + 34 u(MCLin8(2,1,0)) + 34 u(MCLin9(2,1,0))
   - 4 u(MCLinV(1,2,1)) - 7 u(MCLinV(2,1,1)) + 29 u(MCLin3(1,2,1))
   - 4 u(MCLin5(1,2,1)) + 30 u(MCLin8(2,1,1)) + 30 u(MCLin9(2,1,1))
   - 7 u(MCLinV(1,2,2)) - 6 u(MCLinV(2,1,2)) + 36 u(MCLin3(1,2,2))
   - 14 u(MCLin5(1,2,2)) - 13 u(MCLin7(2,1,2)) + 37 u(MCLin9(2,1,2))
   - u(makespan(0)) - 7 u(makespan(1)) - 2 u(makespan(2))
   - u(cut_clique(012,1)) - 41 u(cut_basic_best(02,0))
   - 43 u(cut_two_jobs(0,2,0)) - 13 u(cut_half_best(0,02,0))
   - 14 u(cut_late_job_best(02,0,0,0)) - 16 u(cut_late_job_best(012,0,0,2))
   - 4 u(cut_half_best(1,12,0)) - 4 u(cut_late_job_best(012,0,1,0))
   - 4 u(cut_late_job_best(012,0,1,1)) - 13 u(cut_late_job_best(012,0,1,2))
   - 55 u(cut_basic_plus_epsilon_best(01,0,2))
   - 105 u(cut_basic_best(012,1)) - 49 u(cut_two_jobs(0,2,1))
   - 28 u(cut_basic_plus_epsilon_best(12,1,0))
   - 49 u(cut_basic_plus_epsilon_best(02,1,1))
   - 28 u(cut_basic_plus_epsilon_best(01,1,2))
   - 273 u(cut_basic_best(012,2)) - 147 u(cut_two_jobs(0,1,2))
   - 7 u(cut_half_best(0,012,2))
   - 133 u(cut_basic_plus_epsilon_best(12,2,0))
   - 7 u(cut_late_job_best(012,2,0,0)) - 7 u(cut_late_job_best(012,2,0,1))
   - 7 u(cut_late_job_best(012,2,0,2)) - 13 u(cut_half_best(1,012,2))
   - 133 u(cut_basic_plus_epsilon_best(02,2,1))
   - 13 u(cut_late_job_best(012,2,1,0))
   - 13 u(cut_late_job_best(012,2,1,1))
   - 13 u(cut_late_job_best(012,2,1,2))
   - 147 u(cut_basic_plus_epsilon_best(01,2,2)) - a0 - f0 = -0
Bounds
 u(ord(0,1)) <= 0.99
 u(ord(0,2)) <= 0.99
 u(ord(1,1)) <= 0.99
 u(ord(1,2)) <= 0.99
 u(ord(2,1)) <= 0.99
 u(ord(2,2)) <= 0.99
 -0.99 <= u(triangle2(0,1,0)) <= 0.99
 -0.99 <= u(triangle2(0,1,1)) <= 0.99
 -0.99 <= u(triangle2(0,1,2)) <= 0.99
 -0.99 <= u(triangle2(0,2,0)) <= 0.99
 -0.99 <= u(triangle2(0,2,1)) <= 0.99
 -0.99 <= u(triangle2(0,2,2)) <= 0.99
 -0.99 <= u(triangle2(1,2,0)) <= 0.99
 -0.99 <= u(triangle2(1,2,1)) <= 0.99
 -0.99 <= u(triangle2(1,2,2)) <= 0.99
 -0.99 <= u(MCLinV(0,1,0)) <= 0.99
 -0.99 <= u(MCLinV(1,0,0)) <= 0.99
 u(MCLin1(0,1,0)) <= 0.99
 u(MCLin2(0,1,0)) <= 0.99
 u(MCLin3(0,1,0)) <= 0.99
 u(MCLin4(0,1,0)) <= 0.99
 u(MCLin5(0,1,0)) <= 0.99
 u(MCLin6(1,0,0)) <= 0.99
 u(MCLin7(1,0,0)) <= 0.99
 u(MCLin8(1,0,0)) <= 0.99
 u(MCLin9(1,0,0)) <= 0.99
 u(MCLin10(1,0,0)) <= 0.99
 -0.99 <= u(MCLinV(0,1,1)) <= 0.99
 -0.99 <= u(MCLinV(1,0,1)) <= 0.99
 u(MCLin1(0,1,1)) <= 0.99
 u(MCLin2(0,1,1)) <= 0.99
 u(MCLin3(0,1,1)) <= 0.99
 u(MCLin4(0,1,1)) <= 0.99
 u(MCLin5(0,1,1)) <= 0.99
 u(MCLin6(1,0,1)) <= 0.99
 u(MCLin7(1,0,1)) <= 0.99
 u(MCLin8(1,0,1)) <= 0.99
 u(MCLin9(1,0,1)) <= 0.99
 u(MCLin10(1,0,1)) <= 0.99
 -0.99 <= u(MCLinV(0,1,2)) <= 0.99
 -0.99 <= u(MCLinV(1,0,2)) <= 0.99
 u(MCLin1(0,1,2)) <= 0.99
 u(MCLin2(0,1,2)) <= 0.99
 u(MCLin3(0,1,2)) <= 0.99
 u(MCLin4(0,1,2)) <= 0.99
 u(MCLin5(0,1,2)) <= 0.99
 u(MCLin6(1,0,2)) <= 0.99
 u(MCLin7(1,0,2)) <= 0.99
 u(MCLin8(1,0,2)) <= 0.99
 u(MCLin9(1,0,2)) <= 0.99
 u(MCLin10(1,0,2)) <= 0.99
 -0.99 <= u(MCLinV(0,2,0)) <= 0.99
 -0.99 <= u(MCLinV(2,0,0)) <= 0.99
 u(MCLin1(0,2,0)) <= 0.99
 u(MCLin2(0,2,0)) <= 0.99
 u(MCLin3(0,2,0)) <= 0.99
 u(MCLin4(0,2,0)) <= 0.99
 u(MCLin5(0,2,0)) <= 0.99
 u(MCLin6(2,0,0)) <= 0.99
 u(MCLin7(2,0,0)) <= 0.99
 u(MCLin8(2,0,0)) <= 0.99
 u(MCLin9(2,0,0)) <= 0.99
 u(MCLin10(2,0,0)) <= 0.99
 -0.99 <= u(MCLinV(0,2,1)) <= 0.99
 -0.99 <= u(MCLinV(2,0,1)) <= 0.99
 u(MCLin1(0,2,1)) <= 0.99
 u(MCLin2(0,2,1)) <= 0.99
 u(MCLin3(0,2,1)) <= 0.99
 u(MCLin4(0,2,1)) <= 0.99
 u(MCLin5(0,2,1)) <= 0.99
 u(MCLin6(2,0,1)) <= 0.99
 u(MCLin7(2,0,1)) <= 0.99
 u(MCLin8(2,0,1)) <= 0.99
 u(MCLin9(2,0,1)) <= 0.99
 u(MCLin10(2,0,1)) <= 0.99
 -0.99 <= u(MCLinV(0,2,2)) <= 0.99
 -0.99 <= u(MCLinV(2,0,2)) <= 0.99
 u(MCLin1(0,2,2)) <= 0.99
 u(MCLin2(0,2,2)) <= 0.99
 u(MCLin3(0,2,2)) <= 0.99
 u(MCLin4(0,2,2)) <= 0.99
 u(MCLin5(0,2,2)) <= 0.99
 u(MCLin6(2,0,2)) <= 0.99
 u(MCLin7(2,0,2)) <= 0.99
 u(MCLin8(2,0,2)) <= 0.99
 u(MCLin9(2,0,2)) <= 0.99
 u(MCLin10(2,0,2)) <= 0.99
 -0.99 <= u(MCLinV(1,2,0)) <= 0.99
 -0.99 <= u(MCLinV(2,1,0)) <= 0.99
 u(MCLin1(1,2,0)) <= 0.99
 u(MCLin2(1,2,0)) <= 0.99
 u(MCLin3(1,2,0)) <= 0.99
 u(MCLin4(1,2,0)) <= 0.99
 u(MCLin5(1,2,0)) <= 0.99
 u(MCLin6(2,1,0)) <= 0.99
 u(MCLin7(2,1,0)) <= 0.99
 u(MCLin8(2,1,0)) <= 0.99
 u(MCLin9(2,1,0)) <= 0.99
 u(MCLin10(2,1,0)) <= 0.99
 -0.99 <= u(MCLinV(1,2,1)) <= 0.99
 -0.99 <= u(MCLinV(2,1,1)) <= 0.99
 u(MCLin1(1,2,1)) <= 0.99
 u(MCLin2(1,2,1)) <= 0.99
 u(MCLin3(1,2,1)) <= 0.99
 u(MCLin4(1,2,1)) <= 0.99
 u(MCLin5(1,2,1)) <= 0.99
 u(MCLin6(2,1,1)) <= 0.99
 u(MCLin7(2,1,1)) <= 0.99
 u(MCLin8(2,1,1)) <= 0.99
 u(MCLin9(2,1,1)) <= 0.99
 u(MCLin10(2,1,1)) <= 0.99
 -0.99 <= u(MCLinV(1,2,2)) <= 0.99
 -0.99 <= u(MCLinV(2,1,2)) <= 0.99
 u(MCLin1(1,2,2)) <= 0.99
 u(MCLin2(1,2,2)) <= 0.99
 u(MCLin3(1,2,2)) <= 0.99
 u(MCLin4(1,2,2)) <= 0.99
 u(MCLin5(1,2,2)) <= 0.99
 u(MCLin6(2,1,2)) <= 0.99
 u(MCLin7(2,1,2)) <= 0.99
 u(MCLin8(2,1,2)) <= 0.99
 u(MCLin9(2,1,2)) <= 0.99
 u(MCLin10(2,1,2)) <= 0.99
 u(makespan(0)) <= 0.99
 u(makespan(1)) <= 0.99
 u(makespan(2)) <= 0.99
 u(cut_clique(012,1)) <= 0.99
 u(cut_basic_best(02,0)) <= 0.99
 u(cut_two_jobs(0,2,0)) <= 0.99
 u(cut_half_best(0,02,0)) <= 0.99
 u(cut_late_job_best(02,0,0,0)) <= 0.99
 u(cut_late_job_best(012,0,0,2)) <= 0.99
 u(cut_half_best(1,12,0)) <= 0.99
 u(cut_late_job_best(012,0,1,0)) <= 0.99
 u(cut_late_job_best(012,0,1,1)) <= 0.99
 u(cut_late_job_best(012,0,1,2)) <= 0.99
 u(cut_basic_plus_epsilon_best(01,0,2)) <= 0.99
 u(cut_basic_best(012,1)) <= 0.99
 u(cut_two_jobs(0,2,1)) <= 0.99
 u(cut_half_best(0,012,1)) <= 0.99
 u(cut_basic_plus_epsilon_best(12,1,0)) <= 0.99
 u(cut_late_job_best(012,1,0,0)) <= 0.99
 u(cut_late_job_best(012,1,0,1)) <= 0.99
 u(cut_late_job_best(012,1,0,2)) <= 0.99
 u(cut_half_best(1,012,1)) <= 0.99
 u(cut_basic_plus_epsilon_best(02,1,1)) <= 0.99
 u(cut_late_job_best(012,1,1,0)) <= 0.99
 u(cut_late_job_best(012,1,1,1)) <= 0.99
 u(cut_late_job_best(012,1,1,2)) <= 0.99
 u(cut_basic_plus_epsilon_best(01,1,2)) <= 0.99
 u(cut_basic_best(012,2)) <= 0.99
 u(cut_two_jobs(0,1,2)) <= 0.99
 u(cut_half_best(0,012,2)) <= 0.99
 u(cut_basic_plus_epsilon_best(12,2,0)) <= 0.99
 u(cut_late_job_best(012,2,0,0)) <= 0.99
 u(cut_late_job_best(012,2,0,1)) <= 0.99
 u(cut_late_job_best(012,2,0,2)) <= 0.99
 u(cut_half_best(1,012,2)) <= 0.99
 u(cut_basic_plus_epsilon_best(02,2,1)) <= 0.99
 u(cut_late_job_best(012,2,1,0)) <= 0.99
 u(cut_late_job_best(012,2,1,1)) <= 0.99
 u(cut_late_job_best(012,2,1,2)) <= 0.99
 u(cut_basic_plus_epsilon_best(01,2,2)) <= 0.99
 a(C) >= -100000
 a(x(0,1)) >= -100000
 a(x(0,2)) >= -100000
 a(x(0,0)) >= -100000
 a(x(1,0)) >= -100000
 a(x(1,1)) >= -100000
 a(x(1,2)) >= -100000
 a(x(2,1)) >= -100000
 a(x(2,2)) >= -100000
 a(x(2,0)) >= -100000
 a(y(0,1,1)) >= -100000
 a(y(1,0,1)) >= -100000
 a(y(0,1,2)) >= -100000
 a(y(1,0,2)) >= -100000
 a(y(0,2,1)) >= -100000
 a(y(2,0,1)) >= -100000
 a(y(1,2,0)) >= -100000
 a(y(2,1,0)) >= -100000
 a(y(1,2,1)) >= -100000
 a(y(2,1,1)) >= -100000
 a(v(0,1,0)) >= -100000
 a(v(1,0,0)) >= -100000
 a(u(1,0,0)) >= -100000
 a(v(0,1,1)) >= -100000
 a(v(1,0,1)) >= -100000
 a(u(0,1,1)) >= -100000
 a(u(1,0,1)) >= -100000
 a(v(0,1,2)) >= -100000
 a(v(1,0,2)) >= -100000
 a(u(0,1,2)) >= -100000
 a(u(1,0,2)) >= -100000
 a(v(0,2,0)) >= -100000
 a(v(2,0,0)) >= -100000
 a(u(0,2,0)) >= -100000
 a(u(2,0,0)) >= -100000
 a(v(0,2,1)) >= -100000
 a(v(2,0,1)) >= -100000
 a(u(0,2,1)) >= -100000
 a(u(2,0,1)) >= -100000
 a(v(0,2,2)) >= -100000
 a(v(2,0,2)) >= -100000
 a(u(0,2,2)) >= -100000
 a(u(2,0,2)) >= -100000
 a(v(1,2,0)) >= -100000
 a(v(2,1,0)) >= -100000
 a(u(1,2,0)) >= -100000
 a(u(2,1,0)) >= -100000
 a(v(1,2,1)) >= -100000
 a(v(2,1,1)) >= -100000
 a(u(1,2,1)) >= -100000
 a(u(2,1,1)) >= -100000
 a(v(1,2,2)) >= -100000
 a(v(2,1,2)) >= -100000
 a(u(2,1,2)) >= -100000
 f(C) <= 0.99
 f(x(0,1)) <= 0.99
 f(x(0,2)) <= 0.99
 f(x(0,0)) <= 0.99
 f(x(1,0)) <= 0.99
 f(x(1,1)) <= 0.99
 f(x(1,2)) <= 0.99
 f(x(2,1)) <= 0.99
 f(x(2,2)) <= 0.99
 f(x(2,0)) <= 0.99
 f(y(0,1,1)) <= 0.99
 f(y(1,0,1)) <= 0.99
 f(y(0,1,2)) <= 0.99
 f(y(1,0,2)) <= 0.99
 f(y(0,2,1)) <= 0.99
 f(y(2,0,1)) <= 0.99
 f(y(1,2,0)) <= 0.99
 f(y(2,1,0)) <= 0.99
 f(y(1,2,1)) <= 0.99
 f(y(2,1,1)) <= 0.99
 f(v(0,1,0)) <= 0.99
 f(v(1,0,0)) <= 0.99
 f(u(1,0,0)) <= 0.99
 f(v(0,1,1)) <= 0.99
 f(v(1,0,1)) <= 0.99
 f(u(0,1,1)) <= 0.99
 f(u(1,0,1)) <= 0.99
 f(v(0,1,2)) <= 0.99
 f(v(1,0,2)) <= 0.99
 f(u(0,1,2)) <= 0.99
 f(u(1,0,2)) <= 0.99
 f(v(0,2,0)) <= 0.99
 f(v(2,0,0)) <= 0.99
 f(u(0,2,0)) <= 0.99
 f(u(2,0,0)) <= 0.99
 f(v(0,2,1)) <= 0.99
 f(v(2,0,1)) <= 0.99
 f(u(0,2,1)) <= 0.99
 f(u(2,0,1)) <= 0.99
 f(v(0,2,2)) <= 0.99
 f(v(2,0,2)) <= 0.99
 f(u(0,2,2)) <= 0.99
 f(u(2,0,2)) <= 0.99
 f(v(1,2,0)) <= 0.99
 f(v(2,1,0)) <= 0.99
 f(u(1,2,0)) <= 0.99
 f(u(2,1,0)) <= 0.99
 f(v(1,2,1)) <= 0.99
 f(v(2,1,1)) <= 0.99
 f(u(1,2,1)) <= 0.99
 f(u(2,1,1)) <= 0.99
 f(v(1,2,2)) <= 0.99
 f(v(2,1,2)) <= 0.99
 f(u(2,1,2)) <= 0.99
 a0 >= -100000
 f0 <= 0.99
 obj >= -100000
Generals
 a(C) a(x(0,1)) a(x(0,2)) a(x(0,0)) a(x(1,0)) a(x(1,1)) a(x(1,2)) a(x(2,1))
 a(x(2,2)) a(x(2,0)) a(y(0,1,1)) a(y(1,0,1)) a(y(0,1,2)) a(y(1,0,2))
 a(y(0,2,1)) a(y(2,0,1)) a(y(1,2,0)) a(y(2,1,0)) a(y(1,2,1)) a(y(2,1,1))
 a(v(0,1,0)) a(v(1,0,0)) a(u(1,0,0)) a(v(0,1,1)) a(v(1,0,1)) a(u(0,1,1))
 a(u(1,0,1)) a(v(0,1,2)) a(v(1,0,2)) a(u(0,1,2)) a(u(1,0,2)) a(v(0,2,0))
 a(v(2,0,0)) a(u(0,2,0)) a(u(2,0,0)) a(v(0,2,1)) a(v(2,0,1)) a(u(0,2,1))
 a(u(2,0,1)) a(v(0,2,2)) a(v(2,0,2)) a(u(0,2,2)) a(u(2,0,2)) a(v(1,2,0))
 a(v(2,1,0)) a(u(1,2,0)) a(u(2,1,0)) a(v(1,2,1)) a(v(2,1,1)) a(u(1,2,1))
 a(u(2,1,1)) a(v(1,2,2)) a(v(2,1,2)) a(u(2,1,2)) a0
End
