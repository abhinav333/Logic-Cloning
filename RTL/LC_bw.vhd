----------------------------------------------------------------------------------
-- Company: 
-- Engineer: Abhinav
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
generic (N : integer:=8); --N:width, M: height
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

constant LUT_L_CON: bit_vector :=X"78787878F0F0F0F0";
constant LUT_M_CON: bit_vector :=X"788778870FF00FF0";
constant LUT_N_CON: bit_vector :=X"F0F0F0F0F0F0F0F0";

type TROW_SIG is array (N downto 0) of std_logic_vector(N+5 downto 0);
signal t_row:TROW_SIG:= (others => (others => '0'));

--Include this in test_lut
signal t_row_1:std_logic_vector(N-1 downto 0):=(others => '0');

type CARRY_SIG is array (N-1 downto 0) of std_logic_vector(N+5 downto 0);
signal carry:CARRY_SIG:= (others => (others => '0'));

type LUT_SIG is array (N-1 downto 0) of std_logic_vector(N+5 downto 0);
signal lut_out_0:LUT_SIG:= (others => (others => '0'));
signal lut_out_1:LUT_SIG:= (others => (others => '0'));

signal pp_prod:std_logic_vector(N-1 downto 0):=(others => '0');

constant version:integer:=2;

--function to compute the LUTs in each row
function compute_lut_row1(bit_width : integer; pp_row : integer; version : integer) return integer is
variable comp_result : integer;
begin
  comp_result := 1+pp_row+version;
  if(comp_result > bit_width-1 ) then
    comp_result:=bit_width-1;
  end if;  
  return comp_result;
end function compute_lut_row1;

function assign_a(a : std_logic_vector(N-1 downto 0); pp_row : integer; pp_col:integer; version:integer; bit_width:integer) return std_logic is
variable comp_result : std_logic;
begin
   --Default     
  if(pp_row <= bit_width-version-2) then
   if (pp_col=0 or pp_col=1) then
        comp_result:=a(bit_width-pp_row-version-1);
   else
        comp_result:=a(bit_width-pp_row-version+pp_col-2);
   end if;      
  else
    comp_result:=a(pp_col);--Default
  end if;
  return comp_result;
end function assign_a;

function assign_t_row(t_row_l: std_logic_vector(N+5 downto 0); pp_row : integer; pp_col:integer; version:integer; bit_width:integer) return std_logic is
variable t_row_in : std_logic;
begin
  --t_row_in:=t_row_l(pp_col+1);  --Default
  if(pp_row <= bit_width-version-2) then
        if(pp_col =0) then
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
--Generate LUTS
PP_ROW_GENERATE: for R in 0 to N-1 generate --spans multiplier rows
    FIRST_ROWS: if R < N-1 generate
        PP_COL_GENERATE_F1: for C1 in 0 to compute_lut_row1(N,R,version) generate  --spans multiplier cols
            PP_COL_GEN1:if C1 < compute_lut_row1(N,R,version) generate
            LUTA00X : LUT6_2  generic map(
                                INIT => LUT_L_CON
                            )
                            port map(
                                I5 => '1',
                                I4 => '0',
                                I3 => '0',
                                --I2 => t_row(R)(C1+1),
                                I2 => assign_t_row(t_row(R),R,C1,version,N),
                                --I1 => a(C1),--multiplicand bit
                                I1=>assign_a(a,R,C1,version,N),
                                I0 => b(R),  --multiplier bit
                                O5 => lut_out_1(R)(C1), --PP_in output
                                O6 => lut_out_0(R)(C1) --S out put
                            );
                end generate PP_COL_GEN1;
            PP_COL_GEN2:if C1=compute_lut_row1(N,R,version) and R/=1 generate
            LUTA00X : LUT6_2  generic map(
                                INIT => LUT_M_CON
                            )
                            port map(
                                I5 => '1',
                                I4 => '0',
                                I3 => '0',
                                I2 => t_row_1(R),
                                --I2 => assign_t_row(t_row(R),R,C1,version,N),
                                --I1 => a(C1),--multiplicand bit
                                I1=>assign_a(a,R,C1,version,N),
                                I0 => b(R),  --multiplier bit
                                O5 => lut_out_1(R)(C1), --PP_in output
                                O6 => lut_out_0(R)(C1) --S out put
                            );
                end generate PP_COL_GEN2;
            PP_COL_GEN3:if C1=compute_lut_row1(N,R,version) and R=1 generate
            LUTA00X : LUT6_2  generic map(
                                INIT => LUT_M_CON
                            )
                            port map(
                                I5 => '1',
                                I4 => '0',
                                I3 => '1',
                                --I2 => t_row(R)(C1+1),
                                I2 => t_row_1(R),
                                --I2 => assign_t_row(t_row(R),R,C1,version,N),
                                --I1 => a(C1),--multiplicand bit
                                I1=>assign_a(a,R,C1,version,N),
                                I0 => b(R),  --multiplier bit
                                O5 => lut_out_1(R)(C1), --PP_in output
                                O6 => lut_out_0(R)(C1) --S output
                            );
                end generate PP_COL_GEN3;
           end generate PP_COL_GENERATE_F1;
        end generate FIRST_ROWS ;

        LAST_ROWS: if R=N-1 generate
            PP_COL_GENERATE_F1: for C1 in 0 to N generate  --spans multiplier cols
                P1:if C1 < N-1 generate
                    LUTA00X : LUT6_2  generic map(
                                    INIT => LUT_M_CON
                                )
                                port map(
                                    I5 => '1',
                                    I4 => '0',
                                    I3 => '0',
                                    I2 => t_row(R)(C1+1),
                                    I1 => a(C1),--multiplicand bit
                                    I0 => b(R),  --multiplier bit
                                    O5 => lut_out_1(R)(C1), --PP_in output
                                    O6 => lut_out_0(R)(C1) --S out put
                                );
                    end generate P1;
                P2:if C1 = N-1 generate
                    LUTA00X : LUT6_2  generic map(
                                    INIT => LUT_L_CON
                                )
                                port map(
                                    I5 => '1',
                                    I4 => '0',
                                    I3 => '0',
                                    --I2 => t_row(R)(C1+1),
                                    I2 => t_row_1(R),
                                    I1 => a(C1),--multiplicand bit
                                    I0 => b(R),  --multiplier bit
                                    O5 => lut_out_1(R)(C1), --PP_in output
                                    O6 => lut_out_0(R)(C1) --S output
                                );
                   end generate P2;
                P3:if C1 = N generate
                    LUTA00X : LUT6_2  generic map(
                                    INIT => LUT_N_CON
                                )
                                port map(
                                    I5 => '1',
                                    I4 => '0',
                                    I3 => '0',
                                    I2 => '1',
                                    I1 => '0',--multiplicand bit
                                    I0 => '0',  --multiplier bit
                                    O5 => lut_out_1(R)(C1), --PP_in output
                                    O6 => lut_out_0(R)(C1) --S out put
                                );
      
                end generate P3;
         
               end generate PP_COL_GENERATE_F1;
        end generate LAST_ROWS;   
      end generate PP_ROW_GENERATE ;
        
RR_1:for R in 0 to N-1 generate
FIRST_ROWS:if R<N-1 generate 
CC_1: for CC in 0 to (N/4-1) generate --spans multiplier rows
                            CARRYYX : CARRY4
                            port map(
                                CO => carry(R)(CC*4+4 downto CC*4+1), --only co(3) signal is required
                                --O => result(RC)(CC*4+3 downto CC*4),
                                O => t_row(R+1)(CC*4+3 downto CC*4),
                                CI => carry(R)(CC*4), --attach this to co(3) of previous signal
                                CYINIT => '0',
                                DI => lut_out_1(R)(CC*4+3 downto CC*4) ,
                                S => lut_out_0(R)(CC*4+3 downto CC*4)
                                       
                            );
                      end generate CC_1;
                      t_row_1(R+1)<=carry(R)(compute_lut_row1(N,R,version)+1);
     end generate FIRST_ROWS;
     LAST_ROW:if R=N-1 generate 
CC_1: for CC in 0 to (N/4) generate --spans multiplier rows
                            CARRYYX : CARRY4
                            port map(
                                CO => carry(R)(CC*4+4 downto CC*4+1), --only co(3) signal is required
                                --O => result(RC)(CC*4+3 downto CC*4),
                                O => t_row(R+1)(CC*4+3 downto CC*4),
                                CI => carry(R)(CC*4), --attach this to co(3) of previous signal
                                CYINIT => '0',
                                DI => lut_out_1(R)(CC*4+3 downto CC*4) ,
                                S => lut_out_0(R)(CC*4+3 downto CC*4)
                                          
                            );
                      end generate CC_1;
                      --t_row(R+1)(N)<=carry(R)(N);
     end generate LAST_ROW;                  
                      pp_prod(R)<=t_row(R+1)(0);
    end generate RR_1 ;
p<=t_row(N)(N downto 1) & pp_prod(N-1 downto 0);
end Behavioral;