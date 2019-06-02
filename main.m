warning off;
% params section, insert ur data
W = tf([8.6], [0.0648, 0.54 1]);
Q = [1 0; 0 2];
R = 10;
% end params section
disp("________________SS (input)_______________")
ss(W)
A = ans.A;
B = ans.B;
C = ans.C;
disp("________________P (with rank)_______________")
P = [B A*B]
rank(P)

% intermediate calculation for Riccati equation
% uncomment above for usage
%{
unknK = sym('k',[2 2]);
disp("________________-^K*A_______________")
tmpFirst = -unknK * A;
vpa(tmpFirst)

disp("________________-A^T*^K_______________")
tmpSecond = -A'*unknK;
vpa(tmpSecond)

disp("________________^K*B*r^-1*B^T*^K_______________")
tmpThird = unknK*B*R^-1*B'*unknK;
vpa(tmpThird)

disp("________________Riccati matrix_______________")
eqMatrix = tmpFirst + tmpSecond + tmpThird - Q;
vpa(eqMatrix,3)
%}

disp("________________^K_______________")
G = B*R^-1*B';
circumK = Riccati(A,G,Q)
disp("________________G_______________")
G = R^-1*B'*circumK
disp("________________g (using lqr)_______________")
g = lqr(A,B,Q,R)

% substitution of parameters section
% uncomment above for usage
%{
%}

%backup values
backup.Q = Q;
backup.R = R;
backup.g = g;

g1 = g(1);
g2 = g(2);
q11 = Q(1,1);
q22 = Q(2,2);

% first table
delta = 1;
dMatrixG = [
    g1 g2;
    more(g1, delta) g2;
    less(g1, delta) g2;
    g1 more(g2, delta);
    g1 less(g2, delta);
    more(g1, delta) more(g2, delta);
    less(g1, delta) less(g2, delta);
    more(g1, delta) less(g2, delta);
    less(g1, delta) more(g2, delta);
    ];

sz = [9 3];
varTypes = {'double','double','double'};
varNames = {'g1','g2','J'};
table1 = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
index=1;
disp(' ')
disp("________________STARTED_SUBSTITUTION_FOR_TABLE1_______________")
for row = dMatrixG.'
    try
        g1 = row(1);
        g2 = row(2);
        sim('model_regulator');
        j = J.Data(end);
        r = {g1,g2,j};
        table1(index,:) = r;
        index=index+1;
    catch
        j = J.Data(end);
        r = {g1,g2,j};
        table1(index,:) = r;
        index=index+1;
        fprintf('Inconsistent data in g1 = %s, g2 = %s; skipped.\n', vpa(g1), vpa(g2));
        
    end
    
end
disp(' ')
disp("________________TABLE1_______________")
table1
disp(' ')
disp("________________TABLE1_END_______________")


% restore from backup
g1 = backup.g(1);
g2 = backup.g(2);

% second table

delta = 10;
dMatrixQR = [
    q11 q22 R;
    moreMul(q11, delta) moreMul(q22, delta) R;
    lessDiv(q11, delta) lessDiv(q22, delta) R;
    q11 q22 moreMul(R, delta);
    q11 q22 lessDiv(R, delta);
    moreMul(q11, delta) moreMul(q22, delta) moreMul(R, delta);
    lessDiv(q11, delta) lessDiv(q22, delta) lessDiv(R, delta);
    moreMul(q11, delta) moreMul(q22, delta) lessDiv(R, delta);
    lessDiv(q11, delta) lessDiv(q22, delta) moreMul(R, delta);
    ];

sz = [9 4];
varTypes = {'double','double','double','double'};
varNames = {'g1','g2','J','T'};
table2 = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
index=1;
disp(' ')
disp("________________STARTED_SUBSTITUTION_FOR_TABLE2_______________")
for row = dMatrixQR.'
    try
        % init params
        Q = [row(1) 0; 0 row(2)];
        R = row(3);
        g = lqr(A,B,Q,R);
        
        g1 = g(1);
        g2 = g(2);
        q11 = row(1);
        q22 = row(2);
        
        sim('model_regulator');
        
        j = J.Data(end);
        t = J.Time(end);
        r = {g1,g2,j,t};
        table2(index,:) = r;
        index=index+1;
    catch
        j = J.Data(end);
        t = J.Time(end);
        r = {g1,g2,j,t};
        table2(index,:) = r;
        index=index+1;
        fprintf('Inconsistent data in g1 = %s, g2 = %s; skipped.\n', vpa(g1), vpa(g2));
        
    end
    
end
disp(' ')
disp("________________TABLE2_______________")
table2
disp(' ')
disp("________________TABLE2_END_______________")

% restore from backup
Q = backup.Q;
R = backup.R;
g = backup.g;
g1 = g(1);
g2 = g(2);
q11 = Q(1,1);
q22 = Q(2,2);

