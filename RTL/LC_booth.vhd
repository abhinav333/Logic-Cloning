----------------------------------------------------------------------------------
-- Company: 
-- Engineer:  Abhinav
-- 
-- Create Date: 08/21/2021 07:51:22 AM
-- Design Name: 
-- Module Name: approximate_multiplier - Behavioral
-- Project Name: 
-- Target Devices: 
-- Tool Versions: 
-- Description: 
-- 
-- Dependencies: 
-- 
-- Revision:
-- Revision 0.01 - File Created
-- Additional Comments:
-- 
----------------------------------------------------------------------------------


library IEEE;
use IEEE.STD_LOGIC_1164.ALL;

-- Uncomment the following library declaration if using
-- arithmetic functions with Signed or Unsigned values
--use IEEE.NUMERIC_STD.ALL;

-- Uncomment the following library declaration if instantiating
-- any Xilinx leaf cells in this code.
library UNISIM;
use UNISIM.VComponents.all;

entity approximate_multiplier is
generic (N : integer); --N:width, M: heigth
port(
    a : in std_logic_vector(N-1 downto 0); --width of the multiplier
    b : in std_logic_vector(N-1 downto 0); --heigth of the multiplier
    d : in std_logic_vector(N downto 0); --has always to be zero, only necessary since generic generate
    p : out std_logic_vector(2*N-1 downto 0)
);
end approximate_multiplier;

architecture Behavioral of approximate_multiplier is

component LUT6_2
generic (INIT : bit_vector(63 downto 0) := X"0000000000000000");
port(
    O5 : out std_logic;
    O6 : out std_logic;
    I0 : in std_logic;
    I1 : in std_logic;
    I2 : in std_logic;
    I3 : in std_logic;
    I4 : in std_logic;
    I5 : in std_logic
);
end component;

component CARRY4
port(
    CO : out std_logic_vector(3 downto 0); 
    O : out std_logic_vector(3 downto 0);
    CI : in std_logic;
    CYINIT : in std_logic;
    DI : in std_logic_vector(3 downto 0);
    S : in std_logic_vector(3 downto 0)
);
end component;

--LUT types

constant LUT_A1_CON: bit_vector :=X"C66CC66C0DDF0880";
constant LUT_A_CON: bit_vector :=X"F0C3C3A55A3C3CF0";
constant LUT_A2_CON: bit_vector :=X"3666999333333333";
constant LUT_B_CON: bit_vector :=X"0F3C3C3CC3C3C30F";

type TROW_SIG is array (N/2 downto 0) of std_logic_vector(N+5 downto 0);
signal t_row:TROW_SIG:= (others => (others => '0'));

type TROW_SIG1 is array (N/2 downto 0) of std_logic_vector(1 downto 0);
signal t_row_1:TROW_SIG1:= (others => (others => '0'));

type DI_SIG is array (N/2-1 downto 0) of std_logic_vector(N+5 downto 0);
signal d_row:DI_SIG:= (others => (others => '0'));


type CARRY_SIG is array (N/2-1 downto 0) of std_logic_vector(N+5 downto 0);
signal carry:CARRY_SIG:= (others => (others => '0'));


type LUT_SIG is array (N/2-1 downto 0) of std_logic_vector(N+5 downto 0);
signal lut_out_0:LUT_SIG:= (others => (others => '0'));

type LUT_SIG1 is array (N/2-1 downto 0) of std_logic_vector(1 downto 0);
signal lut_out_1:LUT_SIG1:= (others => (others => '0'));


signal pp_prod:std_logic_vector(N-1 downto 0):=  (others => '0');

signal b_app:std_logic_vector(N downto 0):=  (others => '0');

constant version:integer:=1;

function compute_lut_row1(bit_width : integer; pp_row : integer; version : integer) return integer is
variable comp_result : integer;
begin
  comp_result := 3+2*pp_row+version;
  if(comp_result > bit_width+1 ) then
    comp_result:=bit_width+1; --default
  end if;  
  return comp_result;
end function compute_lut_row1;

--2
function assign_a1(a : std_logic_vector(N-1 downto 0); pp_row : integer; pp_col:integer; version:integer; bit_width:integer) return std_logic is
variable comp_result : std_logic;
begin
  if(pp_row <= (bit_width-version-2)/2) then
   if (pp_col=1 or pp_col=2) then
        comp_result:=a(bit_width-2*pp_row-version);
   else
        comp_result:=a(bit_width-2*pp_row-version+pp_col-2);
   end if;      
  else
    comp_result:=a(pp_col);--Default
  end if;
  return comp_result;
end function assign_a1;

function assign_a2(a : std_logic_vector(N-1 downto 0); pp_row : integer; pp_col:integer; version:integer; bit_width:integer) return std_logic is
variable comp_result : std_logic;
begin
   --Default     
  if(pp_row <= (bit_width-version-2)/2) then
   if (pp_col=1 or pp_col=2) then
        comp_result:=a(bit_width-2*pp_row-version-1);
   else
        comp_result:=a(bit_width-2*pp_row-version+pp_col-3);
   end if;      
  else
    comp_result:=a(pp_col-1);--Default
  end if;
  return comp_result;
end function assign_a2;

--3
function assign_t_row(t_row_l: std_logic_vector(N+5 downto 0); pp_row : integer; pp_col:integer; version:integer; bit_width:integer) return std_logic is
variable t_row_in : std_logic;
begin
   if(pp_row < (bit_width-version)/2) then
         if(pp_col =0 or pp_col=1) then
            t_row_in:=t_row_l(1);
        else    
            t_row_in:=t_row_l(pp_col-1);
        end if;   
    elsif((pp_row = (bit_width-version)/2) and (version rem 2 = 1)) then     
        if(pp_col =0 or pp_col=1) then
            t_row_in:=t_row_l(1);
        else    
            t_row_in:=t_row_l(pp_col);
        end if; 
  else
    t_row_in:=t_row_l(pp_col+1);
  end if;
  return t_row_in;
end function assign_t_row;

begin
--Generate LUTs
t_row_1 (0)(1 downto 0)<="11";
b_app<=b & '0';
PP_ROW_GENERATE: for R in 0 to (N/2-1) generate --spans multiplier rows
    PP_COL_GENERATE_F1: for C1 in 0 to compute_lut_row1(N,R,version) generate  --spans multiplier cols
                FIRST_LUT: if C1=0 generate
                    LUTA00X : LUT6_2  generic map(
                                    INIT => LUT_A1_CON
                                )
                                port map(
                                    I0 => a(C1),
                                    I1 => assign_t_row(t_row(R),R,C1,version,N),
                                    I2 => b_app(2*R),
                                    I3 => b_app(2*R+1),--multiplicand bit
                                    I4 => b_app(2*R+2),  --multiplier bit
                                    I5 => '1',
                                    O5 => lut_out_1(R)(0), --PP_inoutput
                                    O6 => lut_out_0(R)(C1) --S output
                                );
                    pp_prod(2*R)<=lut_out_0(R)(C1);
                    carry(R)(0)<=lut_out_1(R)(0);        
                 end generate FIRST_LUT;
                NEXT_LUT:if C1 > 0 and C1 < compute_lut_row1(N,R,version)-1 generate 
                        LUTA00X : LUT6_2  generic map(
                                INIT => LUT_A_CON
                            )
                            port map(
                                I0 => assign_a2(a,R,C1,version,N),
                                I1 => assign_a1(a,R,C1,version,N),
                                I2 => assign_t_row(t_row(R),R,C1,version,N),
                                I3 => b_app(2*R),
                                I4 => b_app(2*R+1),--multiplicand bit
                                I5 => b_app(2*R+2),  --multiplier bit
                                O6 => lut_out_0(R)(C1) --S out put
                            );
                        d_row(R)(C1-1)<=assign_t_row(t_row(R),R,C1,version,N);
                 end generate NEXT_LUT;
                 NEXT_LUT1:if C1 = compute_lut_row1(N,R,version)-1 generate 
                        LUTA00X : LUT6_2  generic map(
                                INIT => LUT_A2_CON
                            )
                            port map(
                                I0 => a(N-1),
                                I1 => t_row_1(R)(0),
                                I2 => b_app(2*R),
                                I3 => b_app(2*R+1),--multiplicand bit
                                I4 => b_app(2*R+2),  --multiplier bit
                                I5 => '1',
                                O5 => lut_out_1(R)(1), --PP_inoutput
                                O6 => lut_out_0(R)(C1) --S output
                            );
                            d_row(R)(C1-1)<=lut_out_1(R)(1);
                 end generate NEXT_LUT1;
           NEXT_LUT2:if C1 = compute_lut_row1(N,R,version) generate 
                        LUTA00X : LUT6_2  generic map(
                                INIT => LUT_B_CON
                            )
                            port map(
                                I0 => a(N-1),
                                I1 => a(N-1),
                                I2 => t_row_1(R)(1),
                                I3 => b_app(2*R),
                                I4 => b_app(2*R+1),--multiplicand bit
                                I5 => b_app(2*R+2),  --multiplier bit
                                O6 => lut_out_0(R)(C1) --S output
                            );
                            d_row(R)(C1-1)<=t_row_1(R)(1);
                 end generate NEXT_LUT2;
                         
           end generate PP_COL_GENERATE_F1;
end generate PP_ROW_GENERATE ;
        
RR_1:for R in 0 to (N/2-1) generate
CC_1: for CC in 0 to ((compute_lut_row1(N,R,version)-1)/4) generate --spans multiplier rows
                            CARRYYX : CARRY4
                            port map(
                                CO => carry(R)(CC*4+4 downto CC*4+1), --only co(3) signal is required
                                O => t_row(R+1)(CC*4+3 downto CC*4),
                                CI => carry(R)(CC*4), --attach this to co(3) of previous signal
                                CYINIT => '0',
                                DI =>d_row(R)(CC*4+3 downto CC*4),
                                S => lut_out_0(R)(CC*4+4 downto CC*4+1)
                            );
                      end generate CC_1;
                      t_row_1(R+1)(0)<=carry(R)(compute_lut_row1(N,R,version));
                      t_row_1(R+1)(1)<=carry(R)(compute_lut_row1(N,R,version));
                      pp_prod(R*2+1)<=t_row(R+1)(0);
    end generate RR_1;
p<=t_row(N/2)(N downto 1) & pp_prod(N-1 downto 0);
end Behavioral;