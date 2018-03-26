----------------------------------------------------------------------------------
-- Company: 
-- Engineer: 
-- 
-- Create Date: 2018/03/09 14:38:00
-- Design Name: 
-- Module Name: clk_180kHz - Behavioral
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


library IEEE;                       --라이브러리 및 패키지 선언
use IEEE.STD_LOGIC_1164.ALL;

entity clk_180KHz is                --entity : 모듈 이름, 입/출력 port 설정
--  Port ( );
    port(
        clk : in std_logic;
        reset : in std_logic;
        clk_out : out std_logic
    );

end clk_180KHz;

architecture Behavioral of clk_180KHz is            --architecture 
    signal temporal : std_logic;                    --선언부 : 데이터 타입, 신호 및 컴포넌트 등 선언
    signal counter : integer range 0 to 346 :=0;
begin
    frequency_divider : process (reset, clk)     -- 구현부 : 회로 구현
    begin
        if(reset ='1') then                         -- reset 설정
            temporal <= '0';
            counter <= 0;
        elsif rising_edge(clk) then                 -- counter = 346일때 까지 clk의 rising edge를 세겠다!
            if (counter = 346) then
                temporal <= NOT(temporal);          -- counter = 346이면 temporal을 반전
                counter <= 0;
            else
                counter <= counter +1;              -- counter < 346일 때, 1씩 증가
            end if;
        end if;
    end process;
    
    clk_out <= temporal;                            -- temporal 신호를 clk_out으로!
    
end Behavioral;
