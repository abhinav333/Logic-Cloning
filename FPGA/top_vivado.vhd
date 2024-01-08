library IEEE;
use IEEE.STD_LOGIC_1164.ALL;

entity top is
generic (N : integer; M : integer);
port(
topclk : in std_logic;
topa : in std_logic_vector(N-1 downto 0);
topb : in std_logic_vector(M-1 downto 0);
topd : in std_logic_vector(N downto 0);
topp : out std_logic_vector(N+M-1 downto 0)
);
end top;

architecture Behavioral of top is

signal reg_in_a : std_logic_vector(N-1 downto 0);
signal reg_in_b : std_logic_vector(M-1 downto 0);
signal reg_in_d: std_logic_vector(N downto 0);
signal reg_out_p : std_logic_vector(N+M-1 downto 0);

component mult_gen_0 
  PORT (
    CLK : IN STD_LOGIC;
    A : IN STD_LOGIC_VECTOR(3 DOWNTO 0);
    B : IN STD_LOGIC_VECTOR(3 DOWNTO 0);
    P : OUT STD_LOGIC_VECTOR(7 DOWNTO 0)
  );
end component;


--mul : entity work.optimised_accurate
--generic map(
--    N => N,
--    M => M
--)
--port map(
--    a => reg_in_a,
--    b => reg_in_b,
--    d => reg_in_d,
--    p => reg_out_p
--);


begin
vivado_ip: mult_gen_0
port map(
    CLK =>topclk,
    A=>topa,
    B=>topb,
    P=>topp
);
end Behavioral;
