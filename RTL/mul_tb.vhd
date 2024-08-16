----------------------------------------------------------------------------------
-- Company: 
-- Engineer: 
-- 
-- Create Date: 12/06/2021 09:52:07 AM
-- Design Name: 
-- Module Name: mul_tb - Behavioral
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
library UNISIM;
use UNISIM.VComponents.all;
-- Uncomment the following library declaration if using
-- arithmetic functions with Signed or Unsigned values
--use IEEE.NUMERIC_STD.ALL;

-- Uncomment the following library declaration if instantiating
-- any Xilinx leaf cells in this code.
--library UNISIM;
--use UNISIM.VComponents.all;

entity mul_tb is
--  Port ( );
end mul_tb;

architecture Behavioral of mul_tb is
component approximate_multiplier is
generic (N : integer); --N:width, M: heigth
port(
    a : in std_logic_vector(N-1 downto 0); --width of the multiplier
    b : in std_logic_vector(N-1 downto 0); --heigth of the multiplier
    p : out std_logic_vector(2*N-1 downto 0)
);
end component;

constant N : integer := 8;
signal Clk:std_logic :='0';
constant Clk_period : time := 10 ns;
signal a_in : std_logic_vector(N-1 downto 0);
signal b_in : std_logic_vector(N-1 downto 0);
signal p_in : std_logic_vector(2*N-1 downto 0);


begin
uut:approximate_multiplier  
generic map(N => N)
port map(a=>a_in,b=>b_in, p=>p_in );

process
begin
wait for 100 ns;
a_in <= "00010110";
b_in <= "00101100";
wait for 10 ns;
wait;
end process;

end Behavioral;
