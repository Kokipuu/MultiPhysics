# MultiPhysics
1st year master's course, Multiphysics, Group 3
 界面に関する研究


# Virtual Environment
## Build
 ./.fortran/Scripts/activate.ps1

## Exit
deactivate


# Fortran Rule
No distinction is made between upper and lower case letters. 

## The First of Program
program ...

## The end of Program
end program ...

## Ignore implicit type declarations
implicit none


Fortran 主プログラムは，program ... ではじまります。
Fortran は大文字と小文字を区別しません。ここでは原則として全て小文字で書くことにします。
! からはじまる行はコメント行です。コメント行はそれを読む人間のための注釈であり，プログラムの実行時には無視されます。
5行目でまず暗黙の型宣言を無効にします。implicit none と宣言したら，以降はプログラム文内で使用する変数は，全て以下の例のように宣言する必要があります。
次にこれから使う変数の型を real（実数）として宣言します。
7行目では，数値 2.5 を変数 watashi に代入しています。人一倍どころか，人より 2.5 倍の思いがあることを示しています。
8行目では，数値 0.4 を変数 kareshi に代入しています。
9行目では，変数 watashi に代入された数値と，変数 kareshi に代入された数値をかけた答えを，変数 love に代入しています。
10行目では，変数 love の値を出力しています。標準出力（画面）への出力は print に統一します。（print * の * については「書式指定」のところで。）
主プログラムの最後は end program ... で終わります。