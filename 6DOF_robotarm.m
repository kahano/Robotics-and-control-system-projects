%% 6 DOF - Robotic arm
D_1 = simplify(k1);
D_2 = simplify(k2);
D_3 = simplify(k3);
D_matrise = simplify(D_1 + D_2 + D_3);
disp(D_matrise);

%% gravity gk_matrise %%
g1 = diff(p,t1);
g2 = diff(p,t2);
g3 = diff(p,t3);
gk_matrise = sym(zeros(3,1));
gk_matrise(1,1) = g1;
gk_matrise(2,1) = g2;
gk_matrise(3,1) = g3;
disp(gk_matrise);

%% Centrifugal / coriolis (Cijk_matrise)%%

c_111 = 1/2*(diff(D_matrise(1,1), t1) + diff(D_matrise(1,1), t1) - diff(D_matrise(1,1),t1));
c_221 = 1/2*( diff(D_matrise(1,2), t2) + diff(D_matrise(1,2), t2) - diff(D_matrise(2,2),t1));
c_331 = 1/2*( diff(D_matrise(1,3), t3) + diff(D_matrise(1,3), t3) - diff(D_matrise(3,3),t1));
c_121 = 1/2*(diff(D_matrise(1,2), t1) + diff(D_matrise(1,1), t2) - diff(D_matrise(1,2),t1));
c_211 = c_121;
c_131 = 1/2*(diff(D_matrise(3,3), t3) + diff(D_matrise(3,1), t3) - diff(D_matrise(1,3),t3));
c_311 = c_131;
c_321 = 1/2*(diff(D_matrise(1,2), t3) + diff(D_matrise(1,3), t2) - diff(D_matrise(3,2),t1));
c_231 = c_321;

c_112 = 1/2*(diff(D_matrise(2,1), t1) + diff(D_matrise(2,1), t1) - diff(D_matrise(1,1),t2));
c_222 = 1/2*( diff(D_matrise(2,2), t2) + diff(D_matrise(2,2), t2) - diff(D_matrise(2,2),t2));
c_332 = 1/2*( diff(D_matrise(2,3), t3) + diff(D_matrise(2,3), t3) - diff(D_matrise(3,3),t2));
c_122 = 1/2*(diff(D_matrise(2,2), t1) + diff(D_matrise(2,1), t2) - diff(D_matrise(1,2),t2));
c_212 = c_122;
c_132 = 1/2*(diff(D_matrise(2,3), t3) + diff(D_matrise(2,1), t3) - diff(D_matrise(1,3),t2));
c_312 = c_132;
c_322 = 1/2*(diff(D_matrise(2,2), t3) + diff(D_matrise(2,3), t2) - diff(D_matrise(3,2),t2));
c_232 = c_322;

c_113 = 1/2*(diff(D_matrise(3,1), t1) + diff(D_matrise(3,1), t1) - diff(D_matrise(1,1),t3));
c_223 = 1/2*( diff(D_matrise(3,2), t2) + diff(D_matrise(3,2), t2) - diff(D_matrise(2,2),t3));
c_333 = 1/2*( diff(D_matrise(3,3), t3) + diff(D_matrise(3,3), t3) - diff(D_matrise(3,3),t3));
c_123 = 1/2*(diff(D_matrise(3,2), t1) + diff(D_matrise(3,1), t2) - diff(D_matrise(1,2),t3));
c_213 = c_123;
c_133 = 1/2*(diff(D_matrise(3,3), t3) + diff(D_matrise(3,1), t3) - diff(D_matrise(1,3),t3));
c_313 = c_133;
c_323 = 1/2*(diff(D_matrise(3,2), t3) + diff(D_matrise(3,3), t2) - diff(D_matrise(3,2),t3));
c_233 = c_323;

%cijk
c_ijk= c_111+c_221+c_331+c_121+c_211+c_131+c_311+c_131+c_231+c_321+ c_112+c_222+c_332+c_122+c_212+c_132+c_312+c_132+c_232+c_322 + c_113+c_223+c_333+c_123+c_213+c_133+c_313+c_133+c_233+c_323;
C_matrise = sym(zeros(3,3));
%c_[k,j] = sum(c_ijk(q) * qi_dot)
C_matrise(1,1) = [c_111 c_221 c_331]*qi_dot;
C_matrise(1,2) = [c_121 c_221 c_321]*qi_dot;
C_matrise(1,3) = [c_131 c_231 c_331]*qi_dot;
C_matrise(2,1) = [c_112 c_212 c_312]*qi_dot;
C_matrise(2,2) = [c_122 c_222 c_322]*qi_dot;
C_matrise(2,3) = [c_132 c_232 c_332]*qi_dot;
C_matrise(3,1) = [c_113 c_213 c_313]*qi_dot;
C_matrise(3,2) = [c_123 c_223 c_323]*qi_dot;
C_matrise(3,3) = [c_133 c_233 c_333]*qi_dot;
disp(C_matrise)
q_dd = [q1_dd; q2_dd; q3_dd];
